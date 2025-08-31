require "hts"
require "./config"
require "./core/coverage_calculator"
require "./core/cigar"
require "./io/bed_reader"
require "./io/output_manager"
require "./stats/depth_stat"
require "./stats/count_stat"
require "./stats/distribution"

module Depth
  class Runner
    extend Core::CoverageUtils
    extend Stats::Distribution

    def initialize(@config : Configuration)
    end

    def run
      @config.validate!

      # Open BAM/CRAM
      bam = HTS::Bam.open(@config.path, threads: @config.threads)
      begin
        bam.load_index
      rescue ex
        raise "Failed to load index for #{@config.path}: #{ex.message}"
      end

      region = IO.parse_region_str(@config.chrom)
      opts = @config.to_options

      # Create output manager
      output = IO::OutputManager.new(@config.prefix)
      output.create_files(@config.no_per_base, @config.has_regions?)

      begin
        # Get reference sequences from header using hts.cr API
        target_names = bam.header.target_names
        target_lengths = bam.header.target_len
        targets = target_names.map_with_index do |name, i|
          Core::Target.new(name, target_lengths[i].to_i32, i)
        end

        sub_targets = if region
                        targets.select { |t| t.name == region.not_nil!.chrom }
                      else
                        targets
                      end

        # Initialize statistics
        global_dist = Array(Int64).new(512, 0_i64)
        region_dist = Array(Int64).new(512, 0_i64)
        global_stat = Stats::DepthStat.new
        cs = Stats::CountStat(Int32).new(@config.use_median ? 65_536 : 0)

        # Handle window/BED regions
        bed_map : Hash(String, Array(Core::Region))? = nil
        window = @config.window_size
        if bed_path = @config.bed_path
          bed_map = IO.read_bed(bed_path)
        end

        # Create coverage calculator
        calculator = Core::CoverageCalculator.new(bam, opts)

        # Process each target
        sub_targets.each_with_index do |t, idx|
          # Skip if no regions for this chrom and we won't write per-base
          if @config.no_per_base && bed_map && !bed_map.not_nil!.has_key?(t.name)
            next
          end

          coverage = Core::Coverage.new(t.length + 1, 0)
          
          # Determine query region
          query_region = if region && t.name == region.not_nil!.chrom
                           region.not_nil!
                         else
                           Core::Region.new(t.name, 0, 0)
                         end

          tid = calculator.calculate(coverage, query_region)
          next if tid == Core::CoverageResult::ChromNotFound.value

          if tid != Core::CoverageResult::NoData.value
            self.class.cumsum!(coverage)
          end

          # Write per-base intervals
          if output.f_perbase
            if tid == Core::CoverageResult::NoData.value
              output.write_per_base_interval(t.name, 0, t.length, 0)
            else
              self.class.each_depth_interval(coverage) do |(s, e, v)|
                output.write_per_base_interval(t.name, s, e, v)
              end
            end
          end

          # Process regions (window or BED)
          if output.f_regions
            process_regions(t, coverage, tid, window, bed_map, cs, output, region_dist)
          end

          # Process per-chromosome distributions and stats
          if tid != Core::CoverageResult::NoData.value
            self.class.bump_distribution!(global_dist, coverage, 0, coverage.size - 1)
            chrom_stat = Stats::DepthStat.from_array(coverage, 0, coverage.size - 2)
            global_stat = global_stat + chrom_stat
            output.write_summary_line(t.name, chrom_stat)
          end

          # Write distributions
          if f_global = output.f_global
            self.class.write_distribution(f_global.as(::IO), t.name, global_dist)
          end
          if f_region = output.f_region
            self.class.write_distribution(f_region.as(::IO), t.name, region_dist)
            region_dist.fill(0_i64)
          end
          global_dist.fill(0_i64)
        end

      ensure
        output.close_all
      end
    end

    private def process_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32, 
                               window : Int32, bed_map : Hash(String, Array(Core::Region))?, 
                               cs : Stats::CountStat(Int32), output : IO::OutputManager, 
                               region_dist : Array(Int64))
      if window > 0
        process_window_regions(t, coverage, tid, window, cs, output, region_dist)
      else
        process_bed_regions(t, coverage, tid, bed_map, cs, output, region_dist)
      end
    end

    private def process_window_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                     window : Int32, cs : Stats::CountStat(Int32), 
                                     output : IO::OutputManager, region_dist : Array(Int64))
      start = 0
      while start < t.length
        stop = Math.min(start + window, t.length)
        me = 0.0
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median
            cs.clear
            (start...stop).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            len = (stop - start).to_f
            sum = 0_i64
            (start...stop).each { |i| sum += coverage[i] }
            me = (sum.to_f / len)
          end
        end
        output.write_region_line(t.name, start, stop, nil, me)
        if tid != Core::CoverageResult::NoData.value
          idx = [me.to_i, region_dist.size - 1].min
          region_dist[idx] += 1
        end
        start = stop
      end
    end

    private def process_bed_regions(t : Core::Target, coverage : Core::Coverage, tid : Int32,
                                  bed_map : Hash(String, Array(Core::Region))?, 
                                  cs : Stats::CountStat(Int32), output : IO::OutputManager, 
                                  region_dist : Array(Int64))
      regs = bed_map.try(&.[t.name]?) || [] of Core::Region
      regs.each do |r|
        me = 0.0
        if tid != Core::CoverageResult::NoData.value
          if @config.use_median
            cs.clear
            (r.start...Math.min(r.stop, coverage.size)).each { |i| cs.add(coverage[i]) }
            me = cs.median.to_f
          else
            len = (r.stop - r.start).to_f
            sum = 0_i64
            (r.start...Math.min(r.stop, coverage.size)).each { |i| sum += coverage[i] }
            me = len > 0 ? sum.to_f / len : 0.0
          end
        end
        output.write_region_line(t.name, r.start, r.stop, r.name, me)
        if tid != Core::CoverageResult::NoData.value && @config.window_size == 0
          self.class.bump_distribution!(region_dist, coverage, r.start, r.stop)
        end
      end
    end
  end
end

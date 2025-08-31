require "hts"
require "./types"
require "./cigar"

module Depth::Core
  class CoverageCalculator
    extend Cigar
    extend CoverageUtils

    def initialize(@bam : HTS::Bam, @options : Options)
    end

    # Check if record should be filtered out
    private def should_skip_record?(rec) : Bool
      return true if rec.mapq < @options.mapq
      return true if @options.min_frag_len >= 0 && rec.isize.abs < @options.min_frag_len
      return true if rec.isize.abs > @options.max_frag_len
      return true if (rec.flag.value & @options.exclude_flag) != 0
      return true if @options.include_flag != 0 && (rec.flag.value & @options.include_flag) == 0
      
      if @options.read_groups.any?
        rg = rec.aux("RG")
        return true unless rg.is_a?(String) && @options.read_groups.includes?(rg)
      end
      
      false
    end

    # Calculate coverage for a single record
    private def calculate_record_coverage(rec, coverage : Coverage)
      if @options.fast_mode
        start_pos = rec.pos.clamp(0, coverage.size - 1)
        coverage[start_pos] += 1
        endp = rec.endpos
        endp = coverage.size - 1 if endp >= coverage.size
        coverage[endp] -= 1
      elsif @options.fragment_mode
        return if rec.flag.read2? || !rec.flag.proper_pair? || rec.flag.supplementary?
        frag_start = Math.min(rec.pos, rec.mate_pos)
        frag_len = rec.isize.abs
        end_pos = frag_start + frag_len
        end_pos = coverage.size - 1 if end_pos >= coverage.size
        coverage[frag_start] += 1
        coverage[end_pos] -= 1
      else
        self.class.inc_coverage(rec.cigar, rec.pos.to_i32, coverage)
      end
    end

    # Initialize coverage array for a chromosome
    def initialize_coverage_array(coverage : Coverage, chrom_len : Int32)
      coverage.clear
      coverage.concat(Array.new(chrom_len + 1, 0))
    end

    # Process region-specific query
    private def process_region_query(a : Coverage, r : Region, tid : Int32) : {Bool, Int32}
      found = false
      chrom_tid = UNKNOWN_CHROM_TID
      
      q_start = (r.start > 0 ? r.start : 0)
      q_stop = (r.stop > 0 ? r.stop : @bam.header.target_len[tid].to_i32)
      region_str = "#{r.chrom}:#{q_start + 1}-#{q_stop}"
      
      @bam.query(region_str) do |rec|
        unless found
          chrom_len = if tid >= 0
                        @bam.header.target_len[tid].to_i32
                      else
                        chrom_tid = rec.tid
                        @bam.header.target_len[chrom_tid].to_i32
                      end
          initialize_coverage_array(a, chrom_len)
          found = true
        end

        next if should_skip_record?(rec)
        calculate_record_coverage(rec, a)
      end
      
      {found, chrom_tid}
    end

    # Process full BAM scan
    private def process_full_scan(a : Coverage, tid : Int32) : {Bool, Int32}
      found = false
      chrom_tid = UNKNOWN_CHROM_TID
      current_tid = Int32::MIN
      
      @bam.each(copy: false) do |rec|
        # Break if we encounter a different chromosome to prevent array corruption
        if current_tid != Int32::MIN && rec.tid != current_tid
          break
        end
        
        unless found
          chrom_tid = (tid >= 0) ? tid : rec.tid
          current_tid = chrom_tid
          chrom_len = @bam.header.target_len[chrom_tid].to_i32
          initialize_coverage_array(a, chrom_len)
          found = true
        end

        next if should_skip_record?(rec)
        calculate_record_coverage(rec, a)
      end
      
      {found, chrom_tid}
    end

    # Returns: CoverageResult enum values or tid
    def calculate(a : Coverage, r : Region?) : Int32
      # Determine tid
      tid = CoverageResult::ChromNotFound.value
      if r
        tid = @bam.header.get_tid(r.chrom)
      else
        # if no region chrom provided, process everything by iteration
        tid = CoverageResult::AllChroms.value
      end

      return CoverageResult::ChromNotFound.value if tid == CoverageResult::ChromNotFound.value

      found, chrom_tid = if r && tid >= 0
                           process_region_query(a, r, tid)
                         else
                           process_full_scan(a, tid)
                         end

      return CoverageResult::NoData.value unless found
      tid >= 0 ? tid : chrom_tid
    end
  end
end

require "../stats/depth_stat"
require "../config"

module Depth::FileIO
  class OutputManager
    getter f_summary : File?
    getter f_global : File?
    getter f_region : File?
    getter f_perbase : File?
    getter f_regions : File?
    getter f_quantized : File?
    getter f_thresholds : File?
    @header_written = false
    @prefix : String

    def initialize(config : Config)
      @prefix = config.prefix
      label = config.mos_style? ? "mosdepth" : "depth"

      @f_summary = File.open("#{@prefix}.#{label}.summary.txt", "w")
      @f_global = File.open("#{@prefix}.#{label}.global.dist.txt", "w")
      @f_region = config.has_regions? ? File.open("#{@prefix}.#{label}.region.dist.txt", "w") : nil
      @f_perbase = config.no_per_base? ? nil : File.open("#{@prefix}.per-base.bed", "w")
      @f_regions = config.has_regions? ? File.open("#{@prefix}.regions.bed", "w") : nil
      @f_quantized = config.has_quantize? ? File.open("#{@prefix}.quantized.bed", "w") : nil
      @f_thresholds = config.has_thresholds? ? File.open("#{@prefix}.thresholds.bed", "w") : nil
    end

    def write_summary_line(region : String, stat : Depth::Stats::DepthStat)
      return unless @f_summary

      unless @header_written
        @f_summary.not_nil! << ["chrom", "length", "bases", "mean", "min", "max"].join("\t") << '\n'
        @header_written = true
      end

      mean = stat.n_bases > 0 ? stat.sum_depth.to_f / stat.n_bases : 0.0
      pr = (ENV["MOFFDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      mean_str = sprintf("%.#{pr}f", mean)
      minv = stat.min_depth == Int32::MAX ? 0 : stat.min_depth
      # mosdepth uses cumulative depth in the 'bases' column
      @f_summary.not_nil! << [region, stat.n_bases, stat.sum_depth, mean_str, minv, stat.max_depth].join("\t") << '\n'
    end

    # Optionally call at end to add a total line like mosdepth
    def write_summary_total(total : Depth::Stats::DepthStat)
      return unless @f_summary
      mean = total.n_bases > 0 ? total.sum_depth.to_f / total.n_bases : 0.0
      pr = (ENV["MOFFDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      mean_str = sprintf("%.#{pr}f", mean)
      minv = total.min_depth == Int32::MAX ? 0 : total.min_depth
      @f_summary.not_nil! << ["total", total.n_bases, total.sum_depth, mean_str, minv, total.max_depth].join("\t") << "\n"
    end

    def write_per_base_interval(chrom : String, start : Int32, stop : Int32, depth : Int32)
      return unless @f_perbase
      @f_perbase.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << depth << "\n"
    end

    def write_region_stat(chrom : String, start : Int32, stop : Int32, name : String?, value : Float64)
      return unless @f_regions
      pr = (ENV["MOFFDEPTH_PRECISION"]?.try &.to_i?) || (ENV["MOSDEPTH_PRECISION"]?.try &.to_i?) || 2
      val_str = sprintf("%.#{pr}f", value)
      if name
        @f_regions.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << name << "\t" << val_str << "\n"
      else
        @f_regions.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << val_str << "\n"
      end
    end

    def write_quantized_interval(chrom : String, start : Int32, stop : Int32, label : String)
      return unless @f_quantized
      @f_quantized.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << label << "\n"
    end

    def write_thresholds_header(thresholds : Array(Int32))
      return unless @f_thresholds
      @f_thresholds.not_nil! << "#chrom\tstart\tend\tregion"
      thresholds.each { |t| @f_thresholds.not_nil! << "\t#{t}X" }
      @f_thresholds.not_nil! << "\n"
    end

    def write_threshold_counts(chrom : String, start : Int32, stop : Int32, name : String?, counts : Array(Int32))
      return unless @f_thresholds
      @f_thresholds.not_nil! << chrom << "\t" << start << "\t" << stop << "\t"
      @f_thresholds.not_nil! << (name || "unknown")
      counts.each { |c| @f_thresholds.not_nil! << "\t" << c }
      @f_thresholds.not_nil! << "\n"
    end

    def close_all
      [@f_summary, @f_global, @f_region, @f_perbase, @f_regions, @f_quantized, @f_thresholds].each(&.try(&.close))
    end
  end
end

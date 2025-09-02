require "../stats/depth_stat"

module Depth::FileIO
  class OutputManager
    @f_summary : File?
    @f_global : File?
    @f_region : File?
    @f_perbase : File?
    @f_regions : File?
    @f_quantized : File?
    @header_written = false

    def initialize(@prefix : String)
      @f_summary = nil
      @f_global = nil
      @f_region = nil
      @f_perbase = nil
      @f_regions = nil
      @f_quantized = nil
    end

    getter f_summary, f_global, f_region, f_perbase, f_regions, f_quantized

    def create_files(no_per_base : Bool = false, has_regions : Bool = false, has_quantize : Bool = false)
      @f_summary = File.open("#{@prefix}.depth.summary.txt", "w")
      @f_global = File.open("#{@prefix}.depth.global.dist.txt", "w")
      @f_region = has_regions ? File.open("#{@prefix}.depth.region.dist.txt", "w") : nil
      @f_perbase = no_per_base ? nil : File.open("#{@prefix}.per-base.bed", "w")
      @f_regions = has_regions ? File.open("#{@prefix}.regions.bed", "w") : nil
      @f_quantized = has_quantize ? File.open("#{@prefix}.quantized.bed", "w") : nil
    end

    def write_summary_line(region : String, stat : Depth::Stats::DepthStat)
      return unless @f_summary

      unless @header_written
        @f_summary.not_nil! << ["chrom", "length", "sum_depth", "mean", "min", "max"].join("\t") << '\n'
        @header_written = true
      end

      mean = stat.n_bases > 0 ? stat.sum_depth.to_f / stat.n_bases : 0.0
      minv = stat.min_depth == Int32::MAX ? 0 : stat.min_depth
      @f_summary.not_nil! << [region, stat.n_bases, stat.sum_depth, mean, minv, stat.max_depth].join("\t") << '\n'
    end

    def write_per_base_interval(chrom : String, start : Int32, stop : Int32, depth : Int32)
      return unless @f_perbase
      @f_perbase.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << depth << "\n"
    end

    def write_region_stat(chrom : String, start : Int32, stop : Int32, name : String?, value : Float64)
      return unless @f_regions
      if name
        @f_regions.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << name << "\t" << value << "\n"
      else
        @f_regions.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << value << "\n"
      end
    end

    def write_quantized_interval(chrom : String, start : Int32, stop : Int32, label : String)
      return unless @f_quantized
      @f_quantized.not_nil! << chrom << "\t" << start << "\t" << stop << "\t" << label << "\n"
    end

    def close_all
      [@f_summary, @f_global, @f_region, @f_perbase, @f_regions, @f_quantized].each(&.try(&.close))
    end
  end
end

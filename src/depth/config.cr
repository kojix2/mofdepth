require "./core/types"

module Depth
  class Configuration
    property prefix : String = ""
    property path : String = ""
    property threads : Int32 = 0
    property chrom : String = ""
    property by : String = "" # numeric window or BED path
    property no_per_base : Bool = false
    property mapq : Int32 = 0
    property min_frag_len : Int32 = -1
    property max_frag_len : Int32 = -1
    property fast_mode : Bool = false
    property fragment_mode : Bool = false
    property use_median : Bool = false
    property thresholds : Array(Int32) = [] of Int32

    def validate!
      raise ArgumentError.new("BAM/CRAM path is required") if path.empty?
      raise ArgumentError.new("Output prefix is required") if prefix.empty?
      raise ArgumentError.new("Invalid MAPQ threshold") if mapq < 0
      raise ArgumentError.new("Invalid thread count") if threads < 0

      if min_frag_len >= 0 && max_frag_len >= 0 && min_frag_len > max_frag_len
        raise ArgumentError.new("min_frag_len cannot be greater than max_frag_len")
      end

      if fast_mode && fragment_mode
        raise ArgumentError.new("--fast-mode and --fragment-mode cannot be used together")
      end
    end

    def to_options : Core::Options
      Core::Options.new(
        mapq: mapq,
        min_frag_len: min_frag_len,
        max_frag_len: (max_frag_len < 0 ? Int32::MAX : max_frag_len),
        fast_mode: fast_mode,
        fragment_mode: fragment_mode,
      )
    end

    def has_regions? : Bool
      !by.empty?
    end

    def window_size : Int32
      return 0 unless has_regions?
      return 0 unless by.each_char.all?(&.ascii_number?)
      by.to_i
    end

    def bed_path : String?
      return nil unless has_regions?
      return nil if by.each_char.all?(&.ascii_number?)
      by
    end
  end
end

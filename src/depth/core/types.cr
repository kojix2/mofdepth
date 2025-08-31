module Depth::Core
  MAX_COVERAGE      = 400_000
  UNKNOWN_CHROM_TID =    -999

  enum CoverageResult
    ChromNotFound =  -1
    NoData        =  -2
    AllChroms     = -42
  end

  struct Target
    getter name : String
    getter length : Int32
    getter tid : Int32

    def initialize(@name : String, @length : Int32, @tid : Int32)
    end
  end

  struct Region
    getter chrom : String
    getter start : Int32
    getter stop : Int32
    getter name : String?

    def initialize(@chrom : String, @start : Int32, @stop : Int32, @name : String? = nil)
      raise ArgumentError.new("start > stop for #{chrom}:#{start}-#{stop}") if @start > @stop
    end

    def to_s(io)
      if name
        io << chrom << ":" << start + 1 << "-" << stop << "\t" << name
      else
        io << chrom << ":" << start + 1 << "-" << stop
      end
    end
  end

  alias Coverage = Array(Int32)

  record Options,
    mapq : Int32 = 0,
    min_frag_len : Int32 = -1,
    max_frag_len : Int32 = Int32::MAX,
    exclude_flag : UInt16 = 1796_u16,
    include_flag : UInt16 = 0_u16,
    fast_mode : Bool = false,
    fragment_mode : Bool = false,
    read_groups : Array(String) = [] of String
end

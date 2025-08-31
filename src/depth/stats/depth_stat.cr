module Depth::Stats
  struct DepthStat
    property cum_length : Int32 = 0
    property cum_depth : Int64 = 0_i64
    property min_depth : Int32 = Int32::MAX
    property max_depth : Int32 = 0

    def clear
      @cum_length = 0
      @cum_depth = 0
      @min_depth = Int32::MAX
      @max_depth = 0
    end

    def self.from_slice(slice : Slice(Int32))
      s = DepthStat.new
      s.cum_length = slice.size
      slice.each do |v|
        s.cum_depth += v
        s.min_depth = v if v < s.min_depth
        s.max_depth = v if v > s.max_depth
      end
      s
    end

    def self.from_array(array : Array(Int32), start_idx : Int32 = 0, end_idx : Int32 = -1)
      s = DepthStat.new
      end_idx = array.size - 1 if end_idx < 0
      s.cum_length = end_idx - start_idx + 1
      (start_idx..end_idx).each do |i|
        v = array[i]
        s.cum_depth += v
        s.min_depth = v if v < s.min_depth
        s.max_depth = v if v > s.max_depth
      end
      s
    end

    def +(other : DepthStat) : DepthStat
      DepthStat.new.tap do |r|
        r.cum_length = self.cum_length + other.cum_length
        r.cum_depth = self.cum_depth + other.cum_depth
        r.min_depth = {self.min_depth, other.min_depth}.min
        r.max_depth = {self.max_depth, other.max_depth}.max
      end
    end
  end
end

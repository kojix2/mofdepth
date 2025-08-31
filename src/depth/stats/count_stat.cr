module Depth::Stats
  # Histogram-based median for integers
  class CountStat(T)
    getter counts : Array(Int32)
    getter n : Int32 = 0

    def initialize(size : Int32 = 65_536)
      @counts = Array(Int32).new(size, 0)
      @n = 0
    end

    def add(value : Int32)
      @n += 1
      v = value < 0 ? raise ArgumentError.new("negative depth: #{value}") : value
      idx = v < counts.size ? v : counts.size - 1
      counts[idx] += 1
    end

    def median : Int32
      return 0 if n == 0
      stop_n = ((n.to_f * 0.5) + 0.5).to_i
      cum = 0
      counts.each_with_index do |cnt, i|
        cum += cnt
        return i.to_i if cum >= stop_n
      end
      0
    end

    def clear
      return if n == 0
      @n = 0
      counts.fill(0)
    end
  end
end

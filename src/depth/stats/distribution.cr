module Depth::Stats
  module Distribution
    def bump_distribution!(dist : Array(Int64), a : Depth::Core::Coverage, start_i : Int32, stop_i : Int32)
      istop = Math.min(stop_i, a.size)
      i = start_i
      while i < istop
        v = a[i]
        v = Depth::Core::MAX_COVERAGE - 10 if v > Depth::Core::MAX_COVERAGE
        if v >= dist.size
          old = dist.size
          new_size = v + 10
          dist.concat(Array.new(new_size - old, 0_i64))
        end
        dist[v] += 1 if v >= 0
        i += 1
      end
    end

    def write_distribution(io : ::IO, chrom : String, dist : Array(Int64))
      sum = dist.sum
      return if sum < 1
      # cumulative from high → low so line is: depth \t cumulative_fraction
      rev = dist.dup.reverse
      cum = 0.0
      rev.each_with_index do |v, i|
        irev = dist.size - i - 1
        next if irev > 300 && v == 0
        cum += v.to_f / sum
        next if cum < 8e-5
        io << chrom << "\t" << irev << "\t" << cum << '\n'
      end
    end
  end
end

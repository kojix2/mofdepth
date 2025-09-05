module Depth::Core
  module CoverageUtils
    # Utility function for cumulative sum
    def prefix_sum!(a : Coverage)
      sum = 0
      i = 0
      while i < a.size
        sum += a[i]
        a[i] = sum
        i += 1
      end
    end

    # Faster variant: prefix sum only over [start, stop] inclusive, assuming positions < start are already 0.
    # After processing diff-array via this, callers can treat a[0...start] as 0 and a[stop+1...] unchanged.
    def prefix_sum!(a : Coverage, start : Int32, stop : Int32)
      return if a.empty?
      s = start.clamp(0, a.size - 1)
      e = stop.clamp(0, a.size - 1)
      return if e < s

      # compute running sum beginning at s; assumes initial carry is a[s]
      sum = 0
      i = s
      while i <= e
        sum += a[i]
        a[i] = sum
        i += 1
      end
    end

    # Yield intervals [start, stop) with constant depth value
    def each_constant_segment(a : Coverage, stop_at : Int32 = -1, &)
      return if a.size <= 1

      last_depth = Int32::MIN
      last_i = 0
      i = 0
      stop = stop_at <= 0 ? a.size - 1 : stop_at

      while i < stop
        depth = a[i]
        if depth == last_depth
          i += 1
          next
        end

        yield({last_i, i, last_depth}) if last_depth != Int32::MIN
        last_depth = depth
        last_i = i

        break if i + 1 == stop
        i += 1
      end

      if last_i < stop
        yield({last_i, a.size - 1, last_depth})
      elsif last_i != i
        yield({last_i - 1, i, last_depth})
      else
        yield({last_i, i, last_depth})
      end
    end
  end
end

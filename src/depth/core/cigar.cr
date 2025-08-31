require "./types"

module Depth::Core
  module Cigar
    # CIGAR → start/end events on reference
    # Returns array of {pos, +1/-1}; pos is 0-based ref coordinate
    def cigar_start_end_events(cigar, ipos : Int32) : Array(Tuple(Int32, Int32))
      events = [] of Tuple(Int32, Int32)

      pos = ipos
      last_stop = -1
      cigar.each do |op_char, olen|
        # Check if operation consumes reference: M, D, N, =, X
        consumes_ref = ['M', 'D', 'N', '=', 'X'].includes?(op_char)
        next unless consumes_ref

        # Check if operation consumes query: M, I, S, =, X
        consumes_query = ['M', 'I', 'S', '=', 'X'].includes?(op_char)
        if consumes_query
          if pos != last_stop
            events << {pos, 1}
            events << {last_stop, -1} if last_stop >= 0
          end
          last_stop = pos + olen.to_i32
        end
        pos += olen.to_i32
      end
      events << {last_stop, -1} if last_stop >= 0
      events
    end

    def inc_coverage(cigar, ipos : Int32, a : Coverage)
      cigar_start_end_events(cigar, ipos).each do |pos, val|
        next if pos < 0 || a.empty?
        p = pos.clamp(0, a.size - 1)
        a[p] += val
      end
    end
  end

  module CoverageUtils
    # Utility function for cumulative sum
    def cumsum!(a : Coverage)
      sum = 0
      i = 0
      while i < a.size
        sum += a[i]
        a[i] = sum
        i += 1
      end
    end

    # Yield intervals [start, stop) with constant depth value
    def each_depth_interval(a : Coverage, stop_at : Int32 = -1, &)
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

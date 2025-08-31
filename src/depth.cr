require "option_parser"
require "./depth/config"
require "./depth/runner"
require "./depth/version"
require "./depth/errors"

module Depth
  class CLI
    def self.run(args = ARGV)
      config = Configuration.new
      config.prefix = "out"

      parser = OptionParser.parse(args) do |p|
        p.banner = "Usage: depth [options] <prefix> <BAM-or-CRAM>"

        p.on("-t", "--threads THREADS", "BAM decompression threads") { |v| config.threads = v.to_i }
        p.on("-c", "--chrom CHROM", "Restrict to chromosome") { |v| config.chrom = v }
        p.on("-b", "--by BY", "BED file or numeric window size") { |v| config.by = v }
        p.on("-n", "--no-per-base", "Skip per-base output") { config.no_per_base = true }
        p.on("-Q", "--mapq MAPQ", "MAPQ threshold") { |v| config.mapq = v.to_i }
        p.on("-l", "--min-frag-len MIN", "Minimum fragment length") { |v| config.min_frag_len = v.to_i }
        p.on("-u", "--max-frag-len MAX", "Maximum fragment length") { |v| config.max_frag_len = v.to_i }
        p.on("-x", "--fast-mode", "Fast mode") { config.fast_mode = true }
        p.on("-a", "--fragment-mode", "Count full fragment (proper pairs only)") { config.fragment_mode = true }
        p.on("-m", "--use-median", "Use median for region stats instead of mean") { config.use_median = true }
        p.on("-v", "--version", "Show version") { puts Depth::VERSION; exit 0 }
        p.on("-h", "--help", "Show this message") { puts p; exit 0 }

        p.unknown_args do |before, after|
          # last two positional args
          if after.size >= 2
            config.prefix = after[-2]
            config.path = after[-1]
            after.clear
          end
        end
      end

      if ARGV.size < 2
        STDERR.puts "Error: missing arguments: <prefix> <BAM-or-CRAM>"
        STDERR.puts "Use --help for usage information"
        exit 2
      end

      config.prefix = ARGV[0]
      config.path = ARGV[1]

      begin
        runner = Runner.new(config)
        runner.run
      rescue ex : ConfigurationError
        STDERR.puts "Configuration error: #{ex.message}"
        exit 1
      rescue ex : FileNotFoundError
        STDERR.puts "File not found: #{ex.message}"
        exit 1
      rescue ex : IndexError
        STDERR.puts "Index error: #{ex.message}"
        exit 1
      rescue ex : Exception
        STDERR.puts "Error: #{ex.message}"
        exit 1
      end
    end
  end
end

Depth::CLI.run

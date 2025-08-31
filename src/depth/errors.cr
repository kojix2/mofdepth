module Depth
  class Error < Exception; end
  class FileNotFoundError < Error; end
  class InvalidRegionError < Error; end
  class CoverageCalculationError < Error; end
  class ConfigurationError < Error; end
  class IndexError < Error; end
end

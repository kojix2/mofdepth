require "./spec_helper"
require "../src/depth/config"
require "../src/depth/io/output_manager"
require "../src/depth/stats/quantize"

describe "Quantize Integration" do
  describe "OutputManager with quantize" do
    it "creates quantized output file when has_quantize is true" do
      prefix = "test_quantize"
      output = Depth::FileIO::OutputManager.new(prefix)

      # Clean up any existing files
      File.delete("#{prefix}.quantized.bed") if File.exists?("#{prefix}.quantized.bed")

      output.create_files(no_per_base: false, has_regions: false, has_quantize: true)

      # Check that quantized file was created
      output.f_quantized.should_not be_nil

      # Write some test data
      output.write_quantized_interval("chr1", 0, 100, "0:1")
      output.write_quantized_interval("chr1", 100, 200, "1:4")

      output.close_all

      # Verify file exists and has content
      File.exists?("#{prefix}.quantized.bed").should be_true
      content = File.read("#{prefix}.quantized.bed")
      content.should contain("chr1\t0\t100\t0:1")
      content.should contain("chr1\t100\t200\t1:4")

      # Clean up
      File.delete("#{prefix}.quantized.bed")
      File.delete("#{prefix}.depth.summary.txt") if File.exists?("#{prefix}.depth.summary.txt")
      File.delete("#{prefix}.depth.global.dist.txt") if File.exists?("#{prefix}.depth.global.dist.txt")
    end

    it "does not create quantized output file when has_quantize is false" do
      prefix = "test_no_quantize"
      output = Depth::FileIO::OutputManager.new(prefix)

      output.create_files(no_per_base: false, has_regions: false, has_quantize: false)

      # Check that quantized file was not created
      output.f_quantized.should be_nil

      output.close_all

      # Verify file does not exist
      File.exists?("#{prefix}.quantized.bed").should be_false

      # Clean up
      File.delete("#{prefix}.depth.summary.txt") if File.exists?("#{prefix}.depth.summary.txt")
      File.delete("#{prefix}.depth.global.dist.txt") if File.exists?("#{prefix}.depth.global.dist.txt")
    end
  end

  describe "Configuration integration" do
    it "correctly identifies when quantize is enabled" do
      config = Depth::Configuration.new
      config.quantize = "0:1:4:"

      config.has_quantize?.should be_true
      config.quantize_args.should eq([0, 1, 4, Int32::MAX])
    end

    it "correctly identifies when quantize is disabled" do
      config = Depth::Configuration.new
      config.quantize = ""

      config.has_quantize?.should be_false
      config.quantize_args.should eq([] of Int32)
    end
  end

  describe "End-to-end quantize workflow" do
    it "processes quantize workflow correctly" do
      # Test the complete workflow from configuration to output
      config = Depth::Configuration.new
      config.quantize = "0:1:4:"

      # Get quantize args
      quants = config.quantize_args
      quants.should eq([0, 1, 4, Int32::MAX])

      # Create lookup table
      lookup = Depth::Stats::Quantize.make_lookup(quants)
      lookup.should eq(["0:1", "1:4", "4:inf"])

      # Test quantized generation with sample coverage data
      coverage = [0, 0, 1, 1, 1, 4, 4, 0, 0]
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      # Should generate segments based on quantization changes
      segments.size.should be > 0
      segments.each do |start, stop, label|
        start.should be >= 0
        stop.should be > start
        label.should match(/\d+:(inf|\d+)/)
      end
    end

    it "handles edge case with no data correctly" do
      config = Depth::Configuration.new
      config.quantize = "0:1:4:"

      quants = config.quantize_args

      # Test with empty coverage (simulating no data case)
      coverage = [] of Int32
      segments = [] of Tuple(Int32, Int32, String)

      Depth::Stats::Quantize.gen_quantized(quants, coverage) do |segment|
        segments << segment
      end

      # Should handle empty coverage gracefully
      segments.should be_empty
    end

    it "respects environment variables for custom labels" do
      # Set environment variable
      ENV["MOSDEPTH_Q0"] = "LOW"
      ENV["MOSDEPTH_Q1"] = "MEDIUM"
      ENV["MOSDEPTH_Q2"] = "HIGH"

      begin
        quants = [0, 1, 4, Int32::MAX]
        lookup = Depth::Stats::Quantize.make_lookup(quants)

        lookup[0].should eq("LOW")
        lookup[1].should eq("MEDIUM")
        lookup[2].should eq("HIGH")
      ensure
        # Clean up environment variables
        ENV.delete("MOSDEPTH_Q0")
        ENV.delete("MOSDEPTH_Q1")
        ENV.delete("MOSDEPTH_Q2")
      end
    end
  end

  describe "Compatibility with original mosdepth" do
    it "produces same quantize args as original mosdepth test cases" do
      # Test cases from original mosdepth functional-tests.sh

      # Test case: -q 0:1:1000
      config = Depth::Configuration.new
      config.quantize = "0:1:1000"
      args = config.quantize_args
      args.should eq([0, 1, 1000])

      lookup = Depth::Stats::Quantize.make_lookup(args)
      lookup.should eq(["0:1", "1:1000"])
    end

    it "handles single quantize value like original mosdepth" do
      # Test case: -q 60 (single threshold)
      config = Depth::Configuration.new
      config.quantize = "60"
      args = config.quantize_args
      args.should eq([0, 60, Int32::MAX])

      lookup = Depth::Stats::Quantize.make_lookup(args)
      lookup.should eq(["0:60", "60:inf"])
    end
  end
end

require "./spec_helper"
require "file_utils"
require "process"

def run_moffdepth_validation(args : Array(String), prefix : String, test_bam : String, temp_dir : String) : Process::Status
  full_args = args + [prefix, test_bam]
  Process.run("crystal", ["run", "src/depth.cr", "--"] + full_args,
    chdir: Dir.current,
    env: {"PWD" => temp_dir},
    output: Process::Redirect::Close,
    error: Process::Redirect::Close)
end

# Tests to validate moffdepth functionality even in environments where mosdepth doesn't work
describe "moffdepth validation tests" do
  temp_dir = "/tmp/moffdepth_validation"
  mosdepth_dir = "./mosdepth"
  test_bam = "#{mosdepth_dir}/tests/ovl.bam"

  before_each do
    FileUtils.mkdir_p(temp_dir)
  end

  after_each do
    FileUtils.rm_rf(temp_dir) if Dir.exists?(temp_dir)
  end

  describe "basic functionality" do
    it "produces per-base coverage output" do
      status = run_moffdepth_validation([] of String, "#{temp_dir}/test", test_bam, temp_dir)
      status.success?.should be_true

      per_base_file = "#{temp_dir}/test.per-base.bed"
      File.exists?(per_base_file).should be_true
      File.size(per_base_file).should be > 0

      # Check file content structure
      content = File.read(per_base_file)
      lines = content.lines
      lines.size.should be > 0

      # Each line should have 4 columns: chrom, start, end, coverage
      first_line = lines[0].split('\t')
      first_line.size.should eq(4)
    end

    it "produces summary output" do
      status = run_moffdepth_validation([] of String, "#{temp_dir}/summary_test", test_bam, temp_dir)
      status.success?.should be_true

      summary_file = "#{temp_dir}/summary_test.depth.summary.txt"
      File.exists?(summary_file).should be_true
      File.size(summary_file).should be > 0

      content = File.read(summary_file)
      content.should contain("chrom")
      content.should contain("length")
      content.should contain("bases")
    end

    it "produces distribution output" do
      status = run_moffdepth_validation([] of String, "#{temp_dir}/dist_test", test_bam, temp_dir)
      status.success?.should be_true

      dist_file = "#{temp_dir}/dist_test.depth.global.dist.txt"
      File.exists?(dist_file).should be_true
      File.size(dist_file).should be > 0
    end
  end

  describe "flag filtering" do
    it "handles exclude flag filtering" do
      status = run_moffdepth_validation(["-F", "4"], "#{temp_dir}/flag_test", test_bam, temp_dir)
      status.success?.should be_true

      File.exists?("#{temp_dir}/flag_test.per-base.bed").should be_true
    end

    it "handles include flag filtering" do
      status = run_moffdepth_validation(["-i", "2"], "#{temp_dir}/include_test", test_bam, temp_dir)
      status.success?.should be_true

      File.exists?("#{temp_dir}/include_test.per-base.bed").should be_true
    end
  end

  describe "quantization" do
    it "produces quantized output" do
      status = run_moffdepth_validation(["-q", "0:1:1000"], "#{temp_dir}/quant_test", test_bam, temp_dir)
      status.success?.should be_true

      quant_file = "#{temp_dir}/quant_test.quantized.bed"
      File.exists?(quant_file).should be_true
      File.size(quant_file).should be > 0

      # Check quantized format
      content = File.read(quant_file)
      lines = content.lines
      lines.size.should be > 0

      # Each line should have 4 columns: chrom, start, end, quantized_label
      first_line = lines[0].split('\t')
      first_line.size.should eq(4)
      first_line[3].should match(/\d+:(inf|\d+)/)
    end
  end

  describe "region processing" do
    it "handles window-based regions" do
      status = run_moffdepth_validation(["-b", "100"], "#{temp_dir}/window_test", test_bam, temp_dir)
      status.success?.should be_true

      regions_file = "#{temp_dir}/window_test.regions.bed"
      File.exists?(regions_file).should be_true
      File.size(regions_file).should be > 0
    end

    it "handles BED file regions" do
      bed_file = "#{mosdepth_dir}/tests/track.bed"
      next unless File.exists?(bed_file)

      status = run_moffdepth_validation(["-b", bed_file], "#{temp_dir}/bed_test", test_bam, temp_dir)
      status.success?.should be_true

      regions_file = "#{temp_dir}/bed_test.regions.bed"
      File.exists?(regions_file).should be_true
    end
  end

  describe "threshold processing" do
    it "handles threshold calculations" do
      status = run_moffdepth_validation(["-T", "0,1,2", "-b", "100"], "#{temp_dir}/thresh_test", test_bam, temp_dir)
      status.success?.should be_true

      thresh_file = "#{temp_dir}/thresh_test.thresholds.bed"
      File.exists?(thresh_file).should be_true
      File.size(thresh_file).should be > 0

      # Check threshold format
      content = File.read(thresh_file)
      lines = content.lines
      lines.size.should be > 0

      # Should have header + data lines
      header = lines[0]
      header.should contain("0X")
      header.should contain("1X")
      header.should contain("2X")
    end
  end

  describe "advanced modes" do
    it "handles fast mode" do
      status = run_moffdepth_validation(["-x"], "#{temp_dir}/fast_test", test_bam, temp_dir)
      status.success?.should be_true

      File.exists?("#{temp_dir}/fast_test.per-base.bed").should be_true
    end

    it "handles fragment mode" do
      fragment_bam = "#{mosdepth_dir}/tests/full-fragment-pairs.bam"
      next unless File.exists?(fragment_bam)

      full_args = ["-a", "#{temp_dir}/frag_test", fragment_bam]
      status = Process.run("crystal", ["run", "src/depth.cr", "--"] + full_args,
        chdir: Dir.current,
        env: {"PWD" => temp_dir},
        output: Process::Redirect::Close,
        error: Process::Redirect::Close)

      status.success?.should be_true
      File.exists?("#{temp_dir}/frag_test.per-base.bed").should be_true
    end
  end

  describe "error handling" do
    it "handles missing chromosome gracefully" do
      status = run_moffdepth_validation(["-c", "nonexistent"], "#{temp_dir}/missing_test", test_bam, temp_dir)
      # Should fail gracefully
      status.success?.should be_false
    end

    it "handles invalid arguments" do
      status = run_moffdepth_validation(["--invalid-option"], "#{temp_dir}/invalid_test", test_bam, temp_dir)
      # Should fail gracefully
      status.success?.should be_false
    end
  end

  describe "output format validation" do
    it "produces mosdepth-compatible file names" do
      status = run_moffdepth_validation([] of String, "#{temp_dir}/format_test", test_bam, temp_dir)
      status.success?.should be_true

      # Check expected file names match mosdepth format
      File.exists?("#{temp_dir}/format_test.per-base.bed").should be_true
      File.exists?("#{temp_dir}/format_test.depth.summary.txt").should be_true
      File.exists?("#{temp_dir}/format_test.depth.global.dist.txt").should be_true
    end

    it "produces valid BED format output" do
      status = run_moffdepth_validation([] of String, "#{temp_dir}/bed_format_test", test_bam, temp_dir)
      status.success?.should be_true

      per_base_file = "#{temp_dir}/bed_format_test.per-base.bed"
      content = File.read(per_base_file)
      lines = content.lines.reject(&.empty?)

      lines.each do |line|
        fields = line.split('\t')
        fields.size.should eq(4)

        # chrom, start, end, coverage
        chrom = fields[0]
        start_pos = fields[1].to_i
        end_pos = fields[2].to_i
        coverage = fields[3].to_i

        chrom.should_not be_empty
        start_pos.should be >= 0
        end_pos.should be > start_pos
        coverage.should be >= 0
      end
    end
  end
end

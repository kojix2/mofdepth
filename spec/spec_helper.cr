require "spec"
require "process"
require "file_utils"

# Helper to use the compiled binary instead of `crystal run` for speed
module TestBin
  @@built = false

  def self.binary : String
    ENV["MOFFDEPTH_BIN"]? || File.expand_path("./bin/depth")
  end

  def self.ensure_built!
    return if @@built && File.exists?(binary)
    unless File.exists?(binary)
      status = Process.run("shards", ["build", "--release"],
        chdir: Dir.current,
        output: Process::Redirect::Close,
        error: Process::Redirect::Close)
      raise "Failed to build with 'shards build'" unless status.success?
    end
    @@built = true
  end

  def self.run(args : Array(String), prefix : String, bam : String, temp_dir : String) : Process::Status
    ensure_built!
    Process.run(binary, args + [prefix, bam],
      chdir: Dir.current,
      env: {"PWD" => temp_dir},
      output: Process::Redirect::Close,
      error: Process::Redirect::Close)
  end
end

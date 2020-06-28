using PackageCompiler
create_sysimage(:GridapDistributedBenchmark,
  sysimage_path=joinpath(@__DIR__,"..","GridapDistributedBenchmark.so"),
  precompile_execution_file=joinpath(@__DIR__,"..","src","GridapDistributedBenchmark.jl"))

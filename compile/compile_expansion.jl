using PackageCompiler

create_sysimage([:Gridap,:PartitionedArrays,:GridapMHD],
  sysimage_path=joinpath(@__DIR__,"..","GridapMHD.so"),
  precompile_execution_file=joinpath(@__DIR__,"warmup_expansion.jl"))

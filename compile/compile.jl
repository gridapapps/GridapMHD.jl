using PackageCompiler
using Pkg

pkgs = Symbol[]
if VERSION >= v"1.4"
    append!(pkgs, [Symbol(v.name) for v in values(Pkg.dependencies()) if v.is_direct_dep],)
else
    append!(pkgs, [Symbol(name) for name in keys(Pkg.installed())])
end

create_sysimage(pkgs,
  sysimage_path=joinpath(@__DIR__,"..","GridapMHD.so"),
  precompile_execution_file=joinpath(@__DIR__,"warmup.jl"))

  # create_sysimage([:Gridap,:PartitionedArrays,:GridapMHD],
#   sysimage_path=joinpath(@__DIR__,"..","GridapMHD.so"),
#   precompile_execution_file=joinpath(@__DIR__,"warmup.jl"))

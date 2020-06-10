module GridapMHD

using Gridap
using IncompleteLU
using LineSearches: BackTracking

include("InductionlessMHD.jl")

include("LinearSolvers.jl")

using GridapMHD.LinearSolvers

# InductionlessMHD exports
export vprod

# LinearSolvers exports
export GmresSolver
export symbolic_setup
export numerical_setup
export numerical_setup!
export solve!

export writePVD
export compute_u_j_errors

# Benchmarks
include("Hunt.jl")

include("Shercliff.jl")

# Handy functions
_l2(v) = v*v
_h1(v) = v*v + inner(∇(v),∇(v))
_hdiv(v) = v*v + inner((∇*v),(∇*v))

function compute_u_j_errors(uh, jh, u, j, trian, quad)
  eu = uh - u
  ej = jh - j

  eu_l2 = sqrt(sum(integrate(_l2(eu),trian,quad)))
  ej_l2 = sqrt(sum(integrate(_l2(ej),trian,quad)))

  return eu_l2, ej_l2
end

function writePVD(filename,timeSteps)
  rm(filename,force=true,recursive=true)
  mkdir(filename)
  pvdcontent  = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>\n"""
  for t in timeSteps
    pvdcontent *= """    <DataSet timestep=""" * '"'
    pvdcontent *= string(t) * '"' * """ group="" part="0" file=""" * '"'
    pvdcontent *= filename*"""/time_"""*string(t)*""".vtu"/>\n"""
  end
  pvdcontent  *= "  </Collection>\n</VTKFile>"
  f = open(filename * ".pvd", "w")
  write(f,pvdcontent)
  close(f)
end


end # module

module GridapMHD

using Gridap
using GridapODEs
using GridapODEs.TransientFETools
using GridapODEs.ODETools
import GridapODEs.TransientFETools: ∂t

using ForwardDiff
using TimerOutputs
using IncompleteLU
using LineSearches: BackTracking

include("InductionlessMHD.jl")

include("LinearSolvers.jl")

using GridapMHD.LinearSolvers

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

include("ConductiveThinWall.jl")

include("TransientDuctFlow.jl")

# Handy functions
_l2(v, dΩ) = sqrt(sum(∫(v⋅v)*dΩ ))
_h1(v, dΩ) = sqrt(sum(∫(v⋅v + inner(∇(v),∇(v)) )*dΩ ))
_hdiv(v, dΩ) = sqrt(sum(∫(v⋅v + inner((∇⋅v),(∇⋅v)) )*dΩ ))

function compute_u_j_errors(uh, jh, u, j, dΩ)
  eu = uh - u
  ej = jh - j

  eu_l2 = _l2(eu, dΩ)
  ej_l2 = _l2(ej, dΩ)

  return eu_l2, ej_l2
end

function startPVD(filename)
  rm(filename,force=true,recursive=true)
  mkdir(filename)
  pvdheader  = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>\n"""

  io = open(filename * ".pvd", "w")
  write(io,pvdheader)
  close(io)
end

function write_timestep(filename,t,trian;kwargs...)
  pvdcontent = """    <DataSet timestep=""" * '"'
  pvdcontent *= string(t) * '"' * """ group="" part="0" file=""" * '"'
  pvdcontent *= filename*"""/time_"""*string(t)*""".vtu"/>\n"""

  writevtk(trian, filename*"/time_"*string(t)*".vtu";kwargs...)

  io = open(filename * ".pvd", "a")
  write(io,pvdcontent)
  close(io)
end

function closePVD(filename)
  pvdfooter  = "  </Collection>\n</VTKFile>"
  io = open(filename * ".pvd", "a")
  write(io,pvdfooter)
  close(io)
end

end # module

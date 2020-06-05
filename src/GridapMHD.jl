module GridapMHD

using Gridap

include("AnalyticalFunctions.jl")

include("InductionlessMHD.jl")

include("LinearSolvers.jl")

using GridapMHD.AnalyticalSolutions
using GridapMHD.LinearSolvers

# AnalyticalSolutions
export shercliff_u
export shercliff_j
export hunt_u
export hunt_j

# InductionlessMHD
export vprod

# LinearSolvers
export GmresSolver
export symbolic_setup
export numerical_setup
export numerical_setup!
export solve!

export writePVD

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

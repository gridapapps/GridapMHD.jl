module GridapMHD

using Random
using Parameters
using Gridap
using Gridap.Helpers
using Gridap.Algebra
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.Geometry

export ConductingFluid
export MagneticField
export FluidForce
export VelocityBc
export TractionBc
export CurrentBc
export PotentialBc
export FluidForce
export ConductingThinWall

include("Main.jl")

include("Hunt.jl")

end # module

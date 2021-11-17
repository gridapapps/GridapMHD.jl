module GridapMHD

using Parameters
using Gridap

export Material
export VelocityBc
export TractionBc
export InsulatingBc
export ConductingBc
export FluidForce
export WallLaw
export Params

include("Main.jl")

include("Hunt.jl")

end # module

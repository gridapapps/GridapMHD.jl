module PeriodicPoissonTests

include("../src/GridapMHD.jl")
using .GridapMHD
# using GridapMHD

using Gridap
using Test


model = CartesianDiscreteModel((0,1,0,1),(32,32),[1])


end #module

module InductionlessMHDTests

include("../src/GridapMHD.jl")
using .GridapMHD
# using GridapMHD

using Gridap
using Test

uh = GridapMHD.main()

@test 1==1 # check code gets here!


end #module

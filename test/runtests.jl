# using GridapMHD
include("../src/GridapMHD.jl")
using .GridapMHD

# include("../../../.julia/dev/Gridap/src/FESpaces/FESpaces.jl")
# using .FESpaces
using Gridap
using Test

# Dummy test until periodic BC are supported by Gridap

@test 1==1

#@time @testset "PeriodicBC" begin include("PeriodicBCTests.jl") end

#@time @testset "PeriodicPoisson" begin include("PeriodicPoissonTests.jl") end

#@time @testset "Shercliff" begin include("ShercliffTest.jl") end

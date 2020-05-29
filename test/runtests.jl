# using GridapMHD
include("../src/GridapMHD.jl")
using .GridapMHD

# include("../../../.julia/dev/Gridap/src/FESpaces/FESpaces.jl")
# using .FESpaces
using Gridap
using Test


# @time @testset "PeriodicBC" begin include("PeriodicBCTests.jl") end

# @time @testset "PeriodicPoisson" begin include("PeriodicPoissonTests.jl") end

@time @testset "ManufacturedSolMHDTest" begin include("ManufacturedSolMHDTest.jl") end

@time @testset "ConvergenceMHDTest" begin include("ConvergenceMHDTest.jl") end

@time @testset "ConvergenceWithPeriodicBCMHDTest" begin include("ConvergenceWithPeriodicBCMHDTest.jl") end

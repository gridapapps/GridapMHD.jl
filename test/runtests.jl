using GridapMHD

using Gridap
using Test


@time @testset "ManufacturedSolMHDTest" begin include("ManufacturedSolMHDTest.jl") end

@time @testset "ConvergenceMHDTest" begin include("ConvergenceMHDTest.jl") end

@time @testset "ConvergenceWithPeriodicBCMHDTest" begin include("ConvergenceWithPeriodicBCMHDTest.jl") end

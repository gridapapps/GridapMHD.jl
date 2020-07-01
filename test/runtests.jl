using GridapMHD

using Gridap
using Test


@time @testset "ManufacturedSolMHDTest" begin include("ManufacturedSolMHDTest.jl") end

@time @testset "ConvergenceMHDTest" begin include("ConvergenceMHDTest.jl") end

@time @testset "ConvergenceWithPeriodicBCMHDTest" begin include("ConvergenceWithPeriodicBCMHDTest.jl") end

@time @testset "ShercliffTest" begin include("ShercliffTest.jl") end

@time @testset "HuntTest" begin include("HuntTest.jl") end

@time @testset "ThinWallBCTest" begin include("ThinWallBCTest.jl") end

@time @testset "TransientTest" begin include("TransientTest.jl") end

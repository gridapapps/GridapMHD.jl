module GridapMHDSequentialTests

using Test

@time @testset "Hunt" begin include("hunt_tests.jl") end

@time @testset "Expansion" begin include("expansion_tests.jl") end

@time @testset "Cavity" begin include("cavity_tests.jl") end

@time @testset "Transient" begin include("transient_tests.jl") end

end # module

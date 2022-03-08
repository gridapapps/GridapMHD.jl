module GridapMHDTests

using Test

@time @testset "main" begin include("MainTests.jl") end

@time @testset "Hunt" begin include("HuntTests.jl") end

@time @testset "Expansion" begin include("expansion_tests.jl") end

end # module

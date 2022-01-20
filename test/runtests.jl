module GridapMHDTests

using Test

@time @testset "main" begin include("MainTests.jl") end

@time @testset "Hunt" begin include("HuntTests.jl") end

end # module

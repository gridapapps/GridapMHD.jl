module GridapMHDSequentialTests

using Test

@time @testset "Hunt" begin include("hunt_tests.jl") end

@time @testset "Expansion" begin include("expansion_tests.jl") end

@time @testset "Cavity" begin include("cavity_tests.jl") end

@time @testset "Li2019" begin
  include("cavity_li2019_tests.jl")
  include("hunt_li2019_tests.jl")
  include("expansion_li2019_tests.jl")
end

end # module

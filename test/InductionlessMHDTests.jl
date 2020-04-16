module InductionlessMHDTests

include("../src/GridapMHD.jl")
using .GridapMHD
# using GridapMHD

using Gridap
using Test

function run()
    uh = GridapMHD.main()
end

run()

@test 1==1 # check code gets here!


end #module

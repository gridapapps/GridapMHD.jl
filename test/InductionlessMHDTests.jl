module InductionlessMHDTests

include("../src/GridapMHD.jl")
using .GridapMHD
# using GridapMHD

using Gridap
using Test

function run(partition=(4,4,3))
  uh = GridapMHD.main(partition)
end

run()

@test 1==1 # check code gets here!


end #module

module expansion_tests

using GridapMHD

GridapMHD.expansion()
GridapMHD.expansion(np=2,backend=:sequential)
GridapMHD.expansion(np=1,backend=:mpi)


end # module

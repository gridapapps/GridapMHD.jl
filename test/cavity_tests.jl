module cavity_tests

using GridapMHD

GridapMHD.cavity()
GridapMHD.cavity(solver=:block_gmres)
#GridapMHD.cavity(np=(2,2),backend=:sequential)
#GridapMHD.cavity(np=(2,2),backend=:mpi)

end # module

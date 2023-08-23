module cavity_tests

using GridapMHD

GridapMHD.cavity()
#GridapMHD.cavity(np=(2,2),backend=:sequential)
#GridapMHD.cavity(np=(2,2),backend=:mpi)

end # module

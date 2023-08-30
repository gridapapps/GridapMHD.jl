module cavity_tests

using GridapMHD

# Serial, LUSolver
GridapMHD.cavity()

# Serial, GMRES + block preconditioner
GridapMHD.cavity(solver=:block_gmres)

# Sequential, LUSolver
GridapMHD.cavity(np=2,backend=:sequential)

# Sequential, GMRES + block preconditioner
GridapMHD.cavity(np=2,backend=:sequential,solver=:block_gmres)

# Distributed MPI, LUSolver
#GridapMHD.cavity(np=2,backend=:mpi)

end # module

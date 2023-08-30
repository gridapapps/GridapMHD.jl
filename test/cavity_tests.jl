module cavity_tests

using GridapMHD

# Serial, LUSolver
GridapMHD.cavity()

# Serial, GMRES + block preconditioner
GridapMHD.cavity(solver=:block_gmres)

# Sequential, PETSc - SNES + MUMPS
GridapMHD.cavity(np=2,backend=:sequential,solver=:petsc)

# Sequential, GMRES + block preconditioner
GridapMHD.cavity(np=2,backend=:sequential,solver=:block_gmres)

# Distributed MPI, LUSolver
#GridapMHD.cavity(np=2,backend=:mpi)

end # module

module CavityTestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: cavity

# PETSc - SNES + MUMPS
cavity(np=4,nc=(4,4,4),backend=:mpi,solver=:petsc)

end # module
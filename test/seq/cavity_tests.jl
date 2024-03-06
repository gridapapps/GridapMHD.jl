module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR, SparseArrays
using GridapMHD: cavity

# Serial, LUSolver
cavity()

# Sequential, PETSc - SNES + MUMPS
cavity(np=2,backend=:sequential,solver=:petsc)

end # module

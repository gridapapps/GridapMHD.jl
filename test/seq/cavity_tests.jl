module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR, SparseArrays
using GridapMHD: cavity

# Serial, LUSolver
cavity()

end # module

module ExpansionTestsSequential

using GridapMHD: expansion
using GridapPETSc, SparseMatricesCSR

expansion()
expansion(inlet=:shercliff)
expansion(np=2,backend=:sequential)
expansion(np=1,backend=:mpi)

end # module

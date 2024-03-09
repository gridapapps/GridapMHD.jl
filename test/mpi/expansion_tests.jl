module ExpansionTestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: expansion

expansion(np=4,backend=:mpi)

end
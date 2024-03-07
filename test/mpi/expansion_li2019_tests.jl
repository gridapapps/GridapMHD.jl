module ExpansionLi2019TestsMPI

using GridapPETSc
using SparseMatricesCSR
using GridapMHD: expansion

# GMRES + Li2019 with MUMPS for Dj
expansion(np=4,backend=:mpi,solver=:li2019)

end
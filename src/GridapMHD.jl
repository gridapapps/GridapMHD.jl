module GridapMHD

# using Profile
# using FileIO
# using MPI

using Random
using LinearAlgebra
using SparseArrays
using SparseMatricesCSR
using BlockArrays

using FileIO
using BSON

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.MultiField

using GridapDistributed
using GridapPETSc
using GridapPETSc.PETSC
using GridapGmsh

using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector

using PartitionedArrays
using PartitionedArrays: getany

include("Main.jl")

include("ExtraFunctions.jl")

include("Hunt.jl")

include("expansion.jl")

include("cavity.jl")

include("block_solver_li2019.jl")

end # module

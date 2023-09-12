
"""
  This module implements a block-based solver for the MHD equations, 
  based on [(Li,2019)](https://doi.org/10.1137/19M1260372).

  There are three main components: 
    - A non-linear Newton solver, from Gridap.jl 
    - A GMRES linear solver, from GridapSolvers.jl
    - A block-based preconditioner, implemented here and based on [(Li,2019)](https://doi.org/10.1137/19M1260372).
"""
module Li2019

using LinearAlgebra
using BlockArrays

using PartitionedArrays
using Gridap
using GridapDistributed
using GridapPETSc

using GridapSolvers
using GridapSolvers.LinearSolvers
using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector

include("blocks.jl")
include("preconditioner.jl")
include("solver.jl")

export Li2019Solver

end

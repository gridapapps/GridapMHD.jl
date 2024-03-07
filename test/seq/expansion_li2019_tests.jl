
module ExpansionLi2019TestsSequential

using GridapMHD: expansion
using SparseArrays

solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSC{Float64,Int64},
  :vector_type   => Vector{Float64},
  :block_solvers => [:julia,:julia,:cg_jacobi,:cg_jacobi],
  :petsc_options => ""
)

expansion(np=1,backend=:sequential,solver=solver)
expansion(np=2,backend=:sequential,solver=solver)

end # module

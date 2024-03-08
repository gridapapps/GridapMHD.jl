module CavityLi2019TestsSequential

using GridapMHD: cavity
using SparseArrays

solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSC{Float64,Int64},
  :vector_type   => Vector{Float64},
  :block_solvers => [:julia,:julia,:cg_jacobi,:cg_jacobi],
  :petsc_options => ""
)

#cavity(np=1,backend=:sequential,solver=solver)
#cavity(np=2,backend=:sequential,solver=solver)

end # module

module HuntLi2019SequentialTests

using GridapMHD: hunt
using SparseArrays

solver = Dict(
  :solver        => :li2019,
  :matrix_type   => SparseMatrixCSC{0,Float64,Int64},
  :vector_type   => Vector{Float64},
  :block_solvers => [:julia,:julia,:cg_jacobi,:cg_jacobi],
  :petsc_options => ""
)

hunt(
  nc=(4,4),
  np=(1,1),
  backend=:sequential,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=solver,
)

hunt(
  nc=(4,4),
  np=(2,2),
  backend=:sequential,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt",
  solver=solver,
)

end # module
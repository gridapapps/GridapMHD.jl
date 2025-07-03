
using Gridap, GridapDistributed, PartitionedArrays, GridapMHD
using SparseArrays

solver = Dict(
  :solver         => :badia2024,
  :matrix_type    => SparseMatrixCSC{Float64,Int64},
  :vector_type    => Vector{Float64},
  :petsc_options  => "",
  :block_solvers  => [:gmg,:cg_jacobi,:cg_jacobi],
)

hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,10.,0.),
  ζ=10.0,
  order = 2,
  order_j = 3,
  vtk=false,
  title="hunt",
  solver=solver,
  ranks_per_level=fill((1,1,1),2),
  formulation = :mhd,
  rt_scaling = true,
  BL_adapted = false,
)

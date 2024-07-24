module ExpansionBadia2024TestsSequential

using GridapMHD: expansion
using SparseArrays

cw = 0.028
Re = 1.0
Ha = 100.0
N = Ha^2/Re

mesh = Dict(
  :mesher => :p4est_MG,
  :base_mesh => "coarse",
  :num_refs_coarse => 0,
  :ranks_per_level => [1,1],
)
solver = Dict(
  :solver        => :badia2024,
  :matrix_type   => SparseMatrixCSC{Float64,Int64},
  :vector_type   => Vector{Float64},
  :block_solvers => [:gmg,:cg_jacobi,:cg_jacobi],
  :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
)
expansion(np=1,backend=:mpi,mesh=mesh,solver=solver,order=2,ζ=100.0,N=N,Ha=Ha,cw=cw,title="ExpansionGMG",formulation=:mhd)

"""
# ReferenceSolution
expansion(np=1,backend=:mpi,solver=:julia,order=2,ζ=0.0,N=N,Ha=Ha,cw=cw,title="ExpansionRef")

expansion(np=1,backend=:mpi,solver=:julia,order=2,ζ=0.0,N=N,Ha=Ha,cw=cw,title="ExpansionRefBis",mesh=mesh)

meshbis = Dict(
  :mesher => :p4est_SG,
  :base_mesh => "coarse",
  :num_refs => 1
)
expansion(np=1,backend=:mpi,solver=:julia,order=2,ζ=0.0,N=N,Ha=Ha,cw=cw,title="ExpansionRefBis",mesh=meshbis)
"""
end # module
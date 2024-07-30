using Gridap, GridapSolvers, GridapDistributed, PartitionedArrays, GridapPETSc
using BlockArrays, SparseArrays

using Gridap.Algebra, Gridap.FESpaces, Gridap.MultiField
using GridapDistributed: i_am_in
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers

using GridapP4est
using GridapMHD

L = 1.0
u0 = 1.0
ν = 1.0
B = VectorValue(0.0, 10.0, 0.0)
f = VectorValue(0.0, 0.0, 0.0)
B0 = norm(B)
σ = 1.0
ρ = 1.0

Re = u0 * L / ν
Ha = B0 * L * sqrt(σ / (ρ * ν))
N = Ha^2 / Re
f̄ = (L / (ρ * u0^2)) * f
B̄ = (1 / B0) * B
α = 1.0
β = 1.0 / Re
γ = N

# Domain and model
np = 1
ranks = with_mpi() do distribute
  ranks = distribute(LinearIndices((np,)))
end

nc = (6,4,4)
domain = (0, L, 0, L, 0, L)
rank_partition = (np,1,1)
model = CartesianDiscreteModel(ranks,rank_partition,domain, nc)
Ω = Interior(model)

labels = get_face_labeling(model);
Γw = append!(collect(1:4), [9, 10, 13, 14], collect(17:21), collect(23:26))
Γl = append!(collect(5:8), [11, 12, 15, 16, 22])
add_tag_from_tags!(labels, "wall", Γw)
add_tag_from_tags!(labels, "lid", Γl)
add_tag_from_tags!(labels, "insulating", "boundary")
uw = VectorValue(0.0, 0.0, 0.0)
ul = VectorValue(1.0, 0.0, 0.0)
ji = VectorValue(0.0, 0.0, 0.0)

writevtk(model.models.item,"data/cavity_model")

_params = Dict(
  :debug => false,
  :solve => true,
  :res_assemble => false,
  :jac_assemble => false,
  :check_valid => false,
  :model => model,
  :fluid => Dict(
      :domain => model,
      :α => α,
      :β => β,
      :γ => γ,
      :f => f̄,
      :B => B̄,
  ),
  :bcs => Dict(
    :u => Dict(:tags => ["wall", "lid"], :values => [uw, ul]),
    :j => Dict(:tags => "insulating", :values => ji),
  ),
  :solver => :julia,
  :fespaces => Dict(
    :order_u => 2,
    :p_space => :P,
  ),
)

params = GridapMHD.add_default_params(_params);

U, V = GridapMHD._fe_spaces(params)
op = GridapMHD._fe_operator(mfs,U,V,params)

xh = zero(get_trial(op))
solver = GridapMHD._solver(op,params)
xh, cache = solve!(xh,solver,op);

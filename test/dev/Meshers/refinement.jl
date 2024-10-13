
using Gridap, GridapDistributed, PartitionedArrays
using GridapMHD

np = 4
ranks = with_debug() do distribute
  distribute(LinearIndices((np,)))
end

mesh = Dict{Symbol,Any}(
  :mesher    => :gridap_SG,
  :base_mesh => "710",
  :num_refs  => 2,
)

params = Dict{Symbol,Any}()
model = GridapMHD.expansion_mesh(mesh,ranks,params)
smodel = Gridap.Adaptivity.refine(model;refinement_method="simplexify")

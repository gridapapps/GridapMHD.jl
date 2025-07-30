using Gridap
using GridapDistributed
using PartitionedArrays
using GridapMHD

using Gridap.Geometry

parts = (2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(parts),)))
end

model = CartesianDiscreteModel(ranks,parts,(0,1,0,1),(4,4))
ptopo = Geometry.PatchTopology(ReferenceFE{0},model)

Ωp = Geometry.PatchTriangulation(model,ptopo)

Γp = Boundary(Ωp)
Λp = Skeleton(Ωp)

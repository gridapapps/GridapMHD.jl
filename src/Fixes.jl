
#function Base.copy(a::PSparseMatrix)
#  matrix_partition = similar(a.matrix_partition)
#  copy!(matrix_partition, a.matrix_partition)
#  PSparseMatrix(matrix_partition,a.row_partition,a.col_partition)
#end

function Geometry.get_polytopes(model::GridapDistributed.DistributedDiscreteModel)
  polys = map(get_polytopes,local_views(model))
  return PartitionedArrays.getany(polys)
end
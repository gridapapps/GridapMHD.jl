
using DrWatson
using PartitionedArrays
using Gridap, GridapP4est

using Gridap.Arrays, Gridap.Geometry, Gridap.ReferenceFEs

function stretch(x::VectorValue{3})
  α = 1.0
  β = 2.0
  γ = 0.2
  return VectorValue{3}(α⋅x[1],β⋅x[2],γ⋅x[3])
end

function generate_base_mesh()
  vertices = [
    (0,0,0),(1,0,0),(0,1,0),(1,1,0), # Main block, 1st layer
    (0,0,1),(1,0,1),(0,1,1),(1,1,1), # Main block, 2nd layer
    (0,0,2),(1,0,2),(0,1,2),(1,1,2), # Main block, 3rd layer
    (0,0,3),(1,0,3),(0,1,3),(1,1,3), # Main block, 4th layer
    (0,2,1),(1,2,1),                 # Extension,  2nd layer
    (0,2,2),(1,2,2),                 # Extension,  3rd layer
  ]
  vertices = map(x -> stretch(VectorValue{3,Float64}(x...)),vertices)

  c2n_map = Table([
    [ 1, 2, 3, 4, 5, 6, 7, 8],
    [ 5, 6, 7, 8, 9,10,11,12],
    [ 9,10,11,12,13,14,15,16],
    [ 7, 8,17,18,11,12,19,20]
  ])

  cell_type = [1,1,1,1]; polys = [HEX]; orientation = Oriented()
  reffes = map(p->LagrangianRefFE(Float64,p,1),polys)
  grid   = UnstructuredGrid(vertices,c2n_map,reffes,cell_type,orientation)
  model  = UnstructuredDiscreteModel(grid)
  return model
end

function refine_mesh(ranks,model,num_refs)
  ref_model = GridapP4est.with(ranks) do
    ref_model = OctreeDistributedDiscreteModel(ranks,model,num_refs)
    return ref_model
  end
  return ref_model
end

function write_mesh(model)
  trian = Triangulation(model)
  writevtk(trian,datadir("P4est/expansion"))
end

function main(distribute,np)
  ranks  = distribute(LinearIndices((prod(np),)))
  cmodel = generate_base_mesh()
  model  = refine_mesh(ranks,cmodel,3)
  write_mesh(model)
end

with_mpi() do distribute
  main(distribute,(1,1,1))
end

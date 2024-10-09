using Gridap
using GridapGmsh
using PartitionedArrays, GridapDistributed

using Gridap.Geometry, Gridap.Helpers, Gridap.Arrays
using Gridap.Adaptivity

function _get_bounding_box(x::AbstractArray{<:Point})
  pmin,pmax = x[1],x[1]
  for p in x
    pmin = min.(pmin,p)
    pmax = max.(pmax,p)
  end
  pmin,pmax
end

function _get_bounding_box(model::DiscreteModel)
  _get_bounding_box(collect(Geometry.get_node_coordinates(model)))
end

function _get_bounding_box(model::GridapDistributed.DistributedDiscreteModel)
  boxes = map(_get_bounding_box,local_views(model))
  boxes = gather(boxes,destination=:all)
  boxes = map(boxes) do boxes
    if ! isempty(boxes)
      pmin,pmax = boxes[1]
      for (bmin,bmax) in boxes
        pmin = min.(pmin,bmin)
        pmax = max.(pmax,bmax)
      end
      pmin,pmax
    else
      @unreachable
    end
  end
  pmin,pmax = PartitionedArrays.getany(boxes)
  pmin,pmax
end

function top_ha_line(model,Z,n=100)
  δ = 1.e6*eps(Float64)
  pmin,pmax = _get_bounding_box(model)
  xmin,xmax = pmin[1]+δ,pmax[1]-δ
  y1,y2 = pmax[2]-δ, pmax[2]/Z
  line = map( x -> x>0 ? Point(x,y1,δ) : Point(x,y2,δ), range(xmin,xmax,n+1))
  return line
end

ranks = with_mpi() do distribute
  distribute(LinearIndices((4,)))
end

file = "/home/user1/Documents/GridapMHD.jl/meshes/Expansion_710.msh"
model = UnstructuredDiscreteModel(GmshDiscreteModel(ranks,file;has_affine_map=true))
#model = UnstructuredDiscreteModel(GmshDiscreteModel(file;has_affine_map=true))
#model = CartesianDiscreteModel(ranks,(2,2,1),(0,1,0,1,0,1),(10,10,10))

smodel = Adaptivity.refine(model;refinement_method="simplexify")

reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(model,reffe)

u(x) = x[1]
uh = interpolate(u,V)

Ωs = Triangulation(smodel)
uhi = GridapDistributed.DistributedCellField(
  map((uhi,Ωi) -> change_domain(uhi,Ωi,ReferenceDomain()),local_views(uh),local_views(Ωs)),
  Ωs
)

line = top_ha_line(model,4.0,40)
vals = uhi(line)

display(vals)
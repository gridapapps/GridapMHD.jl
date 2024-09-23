module Solid

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.Geometry
using GridapMHD: main
using GridapDistributed
using PartitionedArrays

parts = (2,2,1)
ranks = get_part_ids(SequentialBackend(),parts)

ν=1.0
ρ=1.0
σ=1.0
σ1=0.1
σ2=10.0
B=VectorValue(0.0,50.0,0.0)
f=VectorValue(0.0,0.0,1.0)
u0=1.0
B0=norm(B)

# Reduced quantities
L = 1.0 # Do not rescale by length
Re = u0*L/ν
Ha = B0*L*sqrt(σ/(ρ*ν))
N = Ha^2/Re
f̄ = (L/(ρ*u0^2))*f
B̄ = (1/B0)*B
α = 1.0
β = 1.0/Re
γ = N
σ̄1 = σ1/σ
σ̄2 = σ2/σ

Lt = 2
n = ceil(10*Lt)
const Lf = 1.8
domain = (-Lt/2,Lt/2,-Lt/2,Lt/2,-0.1,0.1)
cells = (n,n,3)
layer(x,a) = sign(x)*abs(x)^(1/a)
xmap((x,y,z)) = VectorValue(layer(x,3),layer(y,3),z)
model = CartesianDiscreteModel(ranks,domain,cells,map=xmap,isperiodic=(false,false,true))
Ω = Interior(model)

# This is not really needed, we know that the maximum tag is the
# same in all parts...
solid_tag =  map_parts(model.models) do model
  labels = get_face_labeling(model)
  cell_entity = labels.d_to_dface_to_entity[end]
  solid = maximum(cell_entity) + 1
end
solid_tag = reduce_all(max,solid_tag,init=0)
solid_1 = get_main_part(solid_tag)
solid_2 = solid_1 + 1
fluid = solid_2 + 1
wall = fluid + 1

function find_tag(xs)
  tol = 1.0e-9
  if all(x->(x[1]>0.5*Lf-tol)||x[1]<-0.5*Lf+tol,xs)
    solid_1
  elseif all(x->(x[2]>0.5*Lf-tol)||x[2]<-0.5*Lf+tol,xs)
    solid_2
  else
    fluid
  end
end

D = num_cell_dims(model)

# Dfaces = get_face_gids(model,D)
#
# cell_entity =  map_parts(model.models,Dfaces.partition) do model,partition
#   labels = get_face_labeling(model)
#   cell_entity = labels.d_to_dface_to_entity[end]
#   grid = get_grid(model)
#   cell_xs = get_cell_coordinates(grid)
#   copyto!(cell_entity,map(find_tag,cell_xs))
#   add_tag!(labels,"solid_1",[solid_1])
#   add_tag!(labels,"solid_2",[solid_2])
#   add_tag!(labels,"solid",[solid_1,solid_2])
#   add_tag!(labels,"fluid",[fluid])
#   # Do we need to make tags consistent? No, we know they already are...
#   topo = get_grid_topology(model)
#   for d in 0:D-1  
#   dface_entity = labels.d_to_dface_to_entity[d+1]
#   dface_cells = get_faces(topo,d,D)
#   cache = array_cache(dface_cells)
#   for dface in 1:length(dface_cells)
#     cells = getindex!(cache,dface_cells,dface)
#     solid_found = false
#     fluid_found = false
#     for cell in cells
#       solid_found = solid_found || cell_entity[cell] in (solid_1,solid_2)
#       fluid_found = fluid_found || cell_entity[cell] == fluid
#     end
#     if solid_found && fluid_found
#       dface_entity[dface] = wall
#     end
#   end
#   end
#   add_tag!(labels,"wall",[wall])
#   insulating_tags = vcat(collect(1:(8+12)),collect((1:4).+(8+12+2)))
#   add_tag_from_tags!(labels,"insulating",insulating_tags)
#   cell_entity[partition.oid_to_lid]
# end

# Separate previous code in two steps: 
# 1) generating tags and 
# 2) generating a cell_field to define a variable conductivity based on tags
# Only step 2 needs to be done if tags are generated with, e.g. gmsh

# 1) Generate tags
# 1.1) Materials: 2 solids and a fluid
# 1.2) Boundaries: non-slip walls (at the interface between a solid and the fluid) 
#      and insulating walls (at exterior faces)
map_parts(model.models) do model
  labels = get_face_labeling(model)
  cell_entity = labels.d_to_dface_to_entity[end]
  grid = get_grid(model)
  cell_xs = get_cell_coordinates(grid)
  copyto!(cell_entity,map(find_tag,cell_xs))
  add_tag!(labels,"solid_1",[solid_1])
  add_tag!(labels,"solid_2",[solid_2])
  add_tag!(labels,"solid",[solid_1,solid_2])
  add_tag!(labels,"fluid",[fluid])
  # Do we need to make tags consistent? No, we know they already are...
  topo = get_grid_topology(model)
  for d in 0:D-1  
  dface_entity = labels.d_to_dface_to_entity[d+1]
  dface_cells = get_faces(topo,d,D)
  cache = array_cache(dface_cells)
  for dface in 1:length(dface_cells)
    cells = getindex!(cache,dface_cells,dface)
    solid_found = false
    fluid_found = false
    for cell in cells
      solid_found = solid_found || cell_entity[cell] in (solid_1,solid_2)
      fluid_found = fluid_found || cell_entity[cell] == fluid
    end
    if solid_found && fluid_found
      dface_entity[dface] = wall
    end
  end
  end
  add_tag!(labels,"wall",[wall])
  insulating_tags = vcat(collect(1:(8+12)),collect((1:4).+(8+12+2)))
  add_tag_from_tags!(labels,"insulating",insulating_tags)
  return nothing
end
writevtk(model,"s_model")

# 2) Get cell entities to define σ
Dfaces = get_face_gids(model,D)
cell_entity =  map_parts(model.models,Dfaces.partition) do model,partition
  labels = get_face_labeling(model)
  cell_entity = labels.d_to_dface_to_entity[end]
  cell_entity[partition.oid_to_lid]
end

function entity_to_σ(entity)
  if entity == solid_1
    σ̄1
  elseif entity == solid_2
    σ̄2
  else
    1.0
  end
end

fields = map_parts(Ω.trians,cell_entity) do trian,cell_entity
  CellField(map(entity_to_σ,cell_entity),trian)
end
σ_Ω = GridapDistributed.DistributedCellField(fields)

# Generate Dict to call GridapMHD
params = Dict(
  :ptimer => PTimer(get_part_ids(sequential,1),verbose=true),
  :debug=>false,
  :res_assemble=>false,
  :jac_assemble=>false,
  :solve=>true,
  :model => model,
  :solid => Dict(:domain=>"solid",:σ=>σ_Ω),
  :fluid => Dict(
    :domain=>"fluid",
    :α=>α,
    :β=>β,
    :γ=>γ,
    :B=>B̄,
    :f=>f̄,
   ),
  :bcs=>Dict(
    :u => Dict(:tags=>"wall"),
    :j => Dict(:tags=>"insulating"),
   ),
)

xh, fullparams = main(params)

ūh,p̄h,j̄h,φ̄h = xh
uh = u0*ūh
ph = (ρ*u0^2)*p̄h
jh = (σ*u0*B0)*j̄h
φh = (u0*B0*L)*φ̄h
dh = ∇⋅jh

writevtk(Ω,"s_Ω",order=2,cellfields=["uh"=>uh,"jh"=>jh,"dh"=>dh,"sigma"=>σ_Ω])


end #module


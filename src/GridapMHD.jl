module GridapMHD

using Gridap
using Gridap.Geometry:get_grid_topology, get_node_coordinates, get_face_nodes

using GridapODEs
using GridapODEs.TransientFETools
using GridapODEs.ODETools
import GridapODEs.TransientFETools: ∂t

using GridapEmbedded
using GridapGmsh

using ForwardDiff
using TimerOutputs
using IncompleteLU
using LineSearches: BackTracking, Static

include("InductionlessMHD.jl")

include("LinearSolvers.jl")

using GridapMHD.LinearSolvers

# LinearSolvers exports
export GmresSolver
export symbolic_setup
export numerical_setup
export numerical_setup!
export solve!

## IO functions exports
export Verbosity
export startPVD
export write_timestep
export closePVD

# Handy functions exports
export compute_u_j_errors
export streching
export cuboid
export add_entity!

# Handy functions
_l2(v, dΩ) = sqrt(sum(∫(v⋅v)*dΩ ))
_h1(v, dΩ) = sqrt(sum(∫(v⋅v + inner(∇(v),∇(v)) )*dΩ ))
_hdiv(v, dΩ) = sqrt(sum(∫(v⋅v + inner((∇⋅v),(∇⋅v)) )*dΩ ))

function compute_u_j_errors(uh, jh, u, j, dΩ)
  eu = uh - u
  ej = jh - j

  eu_l2 = _l2(eu, dΩ)
  ej_l2 = _l2(ej, dΩ)

  return eu_l2, ej_l2
end

# IO functions
struct Verbosity
  items::Dict
  function Verbosity(level::Int = 1)
    items = Dict{String,Bool}()
    if level == 0 # Silent
      items["preprocessSteps"] = false
      items["linearSolver"] = false
      items["nonlinearSolver"] = false
      items["timeStep"] = false
    elseif level == 1 # Standard
      items["preprocessSteps"] = false
      items["linearSolver"] = false
      items["nonlinearSolver"] = true
      items["timeStep"] = true
    elseif level == 2 # Verbose
      items["preprocessSteps"] = false
      items["linearSolver"] = true
      items["nonlinearSolver"] = true
      items["timeStep"] = true
    elseif level == 3 # Debug
      items["preprocessSteps"] = true
      items["linearSolver"] = true
      items["nonlinearSolver"] = true
      items["timeStep"] = true
    else
      @error "Unavailable verbosity level. Choose level 0 to 3, from less to more verbosity"
    end
    new(items)
  end
  function Verbosity(items::Dict)
    new(items)
  end
end

function preprocessStepOutput(v::Verbosity,msg)
  if v.items["preprocessSteps"]
    println(msg)
  end
end

function linearSolverOutput(v::Verbosity)
  return v.items["linearSolver"]
end

function nonlinearSolverOutput(v::Verbosity)
  return v.items["nonlinearSolver"]
end

function timeStepOutput(v::Verbosity,dt,t,it)
  if v.items["timeStep"]
    println("######### TIME STEP $it #########")
    println("Current time: $t")
    println("Time step length: $dt")
    println("#################################")
  end
end

function startPVD(filename)
  rm(filename,force=true,recursive=true)
  mkdir(filename)
  pvdheader  = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>\n"""

  io = open(filename * ".pvd", "w")
  write(io,pvdheader)
  close(io)
end

function write_timestep(filename,t,trian;kwargs...)
  pvdcontent = """    <DataSet timestep=""" * '"'
  pvdcontent *= string(t) * '"' * """ group="" part="0" file=""" * '"'
  pvdcontent *= filename*"""/time_"""*string(t)*""".vtu"/>\n"""

  writevtk(trian, filename*"/time_"*string(t)*".vtu";kwargs...)

  io = open(filename * ".pvd", "a")
  write(io,pvdcontent)
  close(io)
end

function closePVD(filename)
  pvdfooter  = "  </Collection>\n</VTKFile>"
  io = open(filename * ".pvd", "a")
  write(io,pvdfooter)
  close(io)
end


function streching(coord;domain=(0.0,1.0,0.0,1.0,0.0,1.0),factors=(1.0,1.0,1.0),
  dirs=(1,2,3))
ncoord = collect(coord.data)
for (i,dir) in enumerate(dirs)
  ξ0 = domain[i*2-1]
  ξ1 = domain[i*2]
  l =  ξ1 - ξ0
  if ξ0 <= coord[dir] <= ξ1
    ξx = (coord[dir] - ξ0)/l
    if ξx < 1e-12
      ncoord[dir] = ξ0
    else
      ncoord[dir] = ((ξx/sqrt(abs(ξx)))^factors[i])*l + ξ0
    end
  end
end
VectorValue(ncoord)
end


# CSG geometry definitions
function cuboid(;dx=1,dy=1,dz=1,x0=Point(0,0,0),name="cuboid",
  faces=["face$i" for i in 1:6])
  e1 = VectorValue(1,0,0)
  e2 = VectorValue(0,1,0)
  e3 = VectorValue(0,0,1)

  plane1 = plane(x0=x0-0.5*dz*e3,v=-e3,name=faces[1])
  plane2 = plane(x0=x0+0.5*dz*e3,v=+e3,name=faces[2])
  plane3 = plane(x0=x0-0.5*dy*e2,v=-e2,name=faces[3])
  plane4 = plane(x0=x0+0.5*dy*e2,v=+e2,name=faces[4])
  plane5 = plane(x0=x0-0.5*dx*e1,v=-e1,name=faces[5])
  plane6 = plane(x0=x0+0.5*dx*e1,v=+e1,name=faces[6])

  geo12 = intersect(plane1,plane2)
  geo34 = intersect(plane3,plane4)
  geo56 = intersect(plane5,plane6)

  intersect(intersect(geo12,geo34),geo56,name=name)
end

# Function that given a filter set a name to those elements that satisfy the
# conditions of the filter
function add_entity!(model,in,name)
  labels = get_face_labeling(model)
  node_coordinates = get_node_coordinates(model)
  entity = num_entities(labels) + 1
  for d in 0:num_dims(model)-1
    facets = get_face_nodes(model,d)
    for (i,facet) in enumerate(facets)
      coord = sum(node_coordinates[facet])/length(facet)
      if in(coord)
        labels.d_to_dface_to_entity[d+1][i] = entity
      end
    end
  end
  add_tag!(labels,name,[entity])
end


# Benchmarks
include("Hunt.jl")

include("Shercliff.jl")

include("ConductiveThinWall.jl")

include("TransientDuctFlow.jl")

include("DriverInductionlessMHD.jl")

include("TransientDriverInductionlessMHD.jl")

end # module

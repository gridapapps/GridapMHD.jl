module ConvergenceMHDTest

include("../src/GridapMHD.jl")
using .GridapMHD

include("../src/Defaults.jl")

using Gridap
import Gridap: ∇
using LinearAlgebra: tr
using Test

vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

u(x) = VectorValue(x[1]+x[2],-x[1]-x[2],1.0)
p(x) = x[1]+x[2]+x[3]
j(x) = VectorValue(x[1]+x[2]+x[3],x[1]+x[2]+x[3],x[1]+x[2]+x[3])
φ(x) = x[1]+x[2]+x[3]

# ∇u(x) = ∇(u)(x)
# Δu(x) = Δ(u)(x)
# ∇p(x) = ∇(p)(x)
# ∇φ(x) = ∇(φ)(x)
# divu(x) = (∇*u)(x)
# divj(x) = (∇*j)(x)

∇u(x) = TensorValue(1.0,1.0,0.0,-1.0,-1.0,0.0,0.0,0.0,0.0)
Δu(x) = VectorValue(0.0,0.0,0.0)
∇p(x) = VectorValue(1.0,1.0,1.0)
∇φ(x) = VectorValue(1.0,1.0,1.0)
divu(x) = 0.0
divj(x) = 3.0

f_u(x) = - Δu(x) + ∇p(x) #- vprod(j(x),B(x)*B0) + (∇u(x)')*u(x)
f_p(x) = divu(x)
f_j(x) = j(x) + ∇φ(x) - vprod(u(x),B(x)*B0)
f_φ(x) = divj(x)

f(x) = [f_u(x),f_p(x),f_j(x),f_φ(x)]

partition=(3,3,2)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5, sign(x[2])*(abs(x[2])*0.5)^0.5, x[3])

ρ = 1.0
ν = 1.0
σ = 1.0
Re = 1.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 0.0
B0 = 1.0
L = 1.0
N = Ha^2/Re
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
# f(x) = [VectorValue(0.0,0.0,-L^3 * K / Re), 0.0, VectorValue(0.0,0.0,0.0), 0.0]
B(x) = VectorValue(0.0,0.0,0.0)

dirichlet_tags = [append!(collect(1:20),[23,24,25,26]),
                  append!(collect(1:20),[23,24,25,26])]

dirichlet_tags_u = collect(1:26)
dirichlet_tags_j = collect(1:26)


function bc(x)
  return [u(x), j(x)]
end



uh =  GridapMHD.main(
  partition = partition,
  map = map,
  periodic_dir = [],
  domain = (-0.5,0.5,-0.5,0.5,0.0,0.3),
  Δt = 0.0,
  num_time_steps = 2,
  maxit = 3,
  use_dimensionless_formulation = false,
  ν = 1.0,
  ρ = 1.0,
  σ = 1.0,
  L = 1.0,
  Re = 1.0,
  Ha = 0.0,
  f_u = f_u,
  f_p = f_p,
  f_j = f_j,
  f_φ = f_φ,
  B0 = B,
  dirichlet_tags_u = dirichlet_tags_u,
  dirichlet_tags_j = dirichlet_tags_j,
  g_x = bc,
  u0 = u
  )


## Print u0

print_u0 = false
if print_u0
  model = CartesianDiscreteModel((-0.5,0.5,-0.5,0.5,0.0,0.3),partition,[3],map)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags[1])

  Vu = FESpace(
   reffe=:Lagrangian, order=2, valuetype=VectorValue{3,Float64},
   conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

  uh = interpolate_everywhere(Vu,u0)

  trian = Triangulation(model)
  writevtk(trian,"results",nsubcells=2,cellfields=["u"=>uh])
end

end #module

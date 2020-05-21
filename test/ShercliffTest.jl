module ShercliffTest

include("../src/GridapMHD.jl")
using .GridapMHD

include("../src/Defaults.jl")

using Gridap
import Gridap: ∇, Δ
using Test

partition=(15,15,3)
domain = (-0.5,0.5,-0.5,0.5,0.0,0.3)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])

ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 50.0
L = 1.0
N = Ha^2/Re
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1)) # Shercliff's
# K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1)) # Hunt's
f_u(x) = VectorValue(0.0,0.0, -L^3 * K / Re)
B = VectorValue(0.0,Ha,0.0)

dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
dirichlet_tags_j = append!(collect(1:20),[23,24,25,26])

u0(x) = Defaults.shercliff_u(0.5,   # a::Float64,        semi-length of side walls
                             0.5,   # b::Float64,        semi-length of Hartmann walls
                             1.0,   # t_w::Float64,      wall thickness
                             0.0,   # σ_w::Float64,      wall conductivity
                             1.0,   # σ::Float64,        fluid conductivity
                             1.0,   # μ::Float64,        fluid viscosity
                           L^3*K/Re, # grad_pz::Float64,  presure gradient
                             Ha ,   # Ha::Float64,       Hartmann number
                             10 ,   # n::Int,            number of sumands included in Fourier series
                             x)     # x)                 evaluation point

j0(x) = Defaults.shercliff_j(0.5,   # a::Float64,        semi-length of side walls
                             0.5,   # b::Float64,        semi-length of Hartmann walls
                             1.0,   # t_w::Float64,      wall thickness
                             0.0,   # σ_w::Float64,      wall conductivity
                             1.0,   # σ::Float64,        fluid conductivity
                             1.0,   # μ::Float64,        fluid viscosity
                           L^3*K/Re, # grad_pz::Float64,  presure gradient
                             Ha ,   # Ha::Float64,       Hartmann number
                             10 ,   # n::Int,            number of sumands included in Fourier series
                             x)     # x)                 evaluation point
xh, trian, _ = main(
  partition=partition, order=2, domain=domain, periodic_dir=[3],
  dirichlet_tags_u=dirichlet_tags_u, dirichlet_tags_j=dirichlet_tags_j,
  use_dimensionless_formulation=true, f_u=f_u, B=B
  )

uh, ph, jh, φh = xh

writevtk(trian,"results", nsubcells=10,cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
                                                   "u"=>u0, "j"=>j0])

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

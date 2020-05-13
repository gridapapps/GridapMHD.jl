module ShercliffTest

include("../src/GridapMHD.jl")
using .GridapMHD

include("../src/Defaults.jl")

using Gridap
import Gridap: ∇, Δ
using Test

partition=(10,10,3)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])

ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 10.0
L = 1.0
N = Ha^2/Re
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
f(x) = [VectorValue(0.0,0.0,-L^3 * K / Re), 0.0, VectorValue(0.0,0.0,0.0), 0.0]
B(x) = VectorValue(0.0,Ha,0.0)

dirichlet_tags = [append!(collect(1:20),[23,24,25,26]),
                  append!(collect(1:20),[23,24,25,26])]

u0(x) = Defaults.shercliff_solution(0.5,   # a::Float64,        semi-length of side walls
                                    0.5,   # b::Float64,        semi-length of Hartmann walls
                                    1.0,   # t_w::Float64,      wall thickness
                                    0.0,   # σ_w::Float64,      wall conductivity
                                    1.0,   # σ::Float64,        fluid conductivity
                                    1.0,   # μ::Float64,        fluid viscosity
                                -L^3*K/Re, # grad_pz::Float64,  presure gradient
                                    Ha ,   # Ha::Float64,       Hartmann number
                                    10 ,   # n::Int,            number of sumands included in Fourier series
                                    x)[1]  # x)                 evaluation point

uh =  main(partition           = partition,
           map                 = map,
           periodic_dir        = [3],
           domain              = (-0.5,0.5,-0.5,0.5,0.0,0.3),
           Δt                  = 1e-2,
           num_time_steps      = 2,
           maxit               = 5,
  use_dimensionless_formulation= true,
           ν                   = 1.0,
           ρ                   = 1.0,
           σ                   = 1.0,
           L                   = 1.0,
           Re                  = 1.0,
           Ha                  = 1.0,
           f                   = f,
           B0                  = B,
           dirichlet_tags      = dirichlet_tags,
           u0                  = u0
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

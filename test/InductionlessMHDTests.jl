module InductionlessMHDTests

include("../src/GridapMHD.jl")
using .GridapMHD
# using GridapMHD
using Gridap
import Gridap: ∇, Δ
using Test

partition=(11,11,3)
uh = GridapMHD.main(partition=partition, Δt=5.0e-3,nt=20,
              use_dimensionless_formulation=true)

# ρ = 1.0
# ν = 1.0
# σ = 1.0
# Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
# Ha = 10.0
# N = σ*L*(Ha^2)/(ρ*Re)
# K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
# f(x) = VectorValue(0.0,0.0,L^3 * K / Re)
# B(x) = VectorValue(0.0,Ha,0.0)

# domain = (-0.5,0.5,-0.5,0.5,0,0.1)
# order = 2
# map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])
# model = CartesianDiscreteModel(domain,partition,[3],map)
# trian = Triangulation(model)
#
# labels = get_face_labeling(model)
# add_tag_from_tags!(labels,"dirichlet",append!(collect(1:20),[23,24,25,26]))
#
# Vu = FESpace(
#  reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
#  conformity=:H1, model=model, dirichlet_tags="dirichlet")
#
# Vj = FESpace(
#  reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
#  conformity=:Hdiv, model=model, dirichlet_tags="dirichlet")
#
# Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
# Ha = 10.0
# K = Ha / (1-0.825*Ha^(1/2)-Ha^(-1))
# grad_pz = K / Re
# u(x) = analytical_solution(0.5,   # a::Float64,        semi-length of side walls
#                       0.5,   # b::Float64,        semi-length of Hartmann walls
#                       1.0,   # t_w::Float64,      wall thickness
#                       0.0,   # σ_w::Float64,      wall conductivity
#                       1.0,   # σ::Float64,        fluid conductivity
#                       1.0,   # μ::Float64,        fluid viscosity
#                     grad_pz, # grad_pz::Float64,  presure gradient
#                       Ha ,   # Ha::Float64,       Hartmann number
#                       10 ,   # n::Int,            number of sumands included in Fourier series
#                       x)[1]  # x)                 evaluation point
# ua(x) = u(x).array
#
# j(x) = analytical_solution(0.5,   # a::Float64,        semi-length of side walls
#                            0.5,   # b::Float64,        semi-length of Hartmann walls
#                            1.0,   # t_w::Float64,      wall thickness
#                            0.0,   # σ_w::Float64,      wall conductivity
#                            1.0,   # σ::Float64,        fluid conductivity
#                            1.0,   # μ::Float64,        fluid viscosity
#                          grad_pz, # grad_pz::Float64,  presure gradient
#                            Ha ,   # Ha::Float64,       Hartmann number
#                            10 ,   # n::Int,            number of sumands included in Fourier series
#                            x)[2]  # x)                 evaluation point
#
# vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])
#
# ∇p(x) = VectorValue(0,0,grad_pz)
# f(x) =  VectorValue(0,0,grad_pz)
# ρ = 1.0
# ν = 1.0
# B(x) = VectorValue(0.0,10.0,0.0)
# ∇u(x) = ∇(u)(x)
# Δu(x) = ∇*(∇(u)(x))
#
# x = get_physical_coordinate(trian)
# jh = interpolate_everywhere(Vj,j)
# uh = interpolate_everywhere(Vu,u)
# ∇uh = gradient(uh)
# jxB = interpolate_everywhere(Vu,x->vprod(j(x),B(x)))
# Δuh = interpolate_everywhere(Vu,Δu)
#
# res_u = uh*(∇uh) - ν*Δuh + ∇p - 1/ρ * jxB - f
#
# res = integrate(res_u,trian,quad)
#
#
# writevtk(trian,"results",cellfields=["u"=>uh])


end #module

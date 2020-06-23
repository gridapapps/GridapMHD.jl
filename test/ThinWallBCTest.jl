module ThinWallBCTest

using Gridap
using Test
using LineSearches: BackTracking

include("../src/GridapMHD.jl")
using .GridapMHD
using .GridapMHD:analytical_shercliff_u
using .GridapMHD:analytical_shercliff_j
include("../src/InductionlessMHD.jl")
using .InductionlessMHD

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0
wall_thickness = 0.10
wall_conductivity = 1.0
c_w = wall_conductivity * wall_thickness/σ

U0 = 10.0
B0 = 10.0
L = 1.0
Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))

nx = 15; ny = 15

N = Ha^2/Re
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3

f_u(x) = VectorValue(0.0,0.0, -∂p∂z) * L/U0^2
g_u = VectorValue(0.0,0.0,0.0)
g_j = VectorValue(0.0,0.0,0.0)
B = VectorValue(0.0,Ha,0.0)/B0

# Discretizatoin
order = 2
domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])*2/sqrt(2)


dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
dirichlet_tags_j = append!(collect(1:20),[25,26])
neumann_j = [23,24]

partition=(nx,ny,3)
model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false,true),map=map)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)
add_tag_from_tags!(labels,"neumann_j",neumann_j)

Vu = FESpace(
  reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
  conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

Vp = FESpace(
  reffe=:PLagrangian, order=order-1, valuetype=Float64,
  conformity=:L2, model=model, constraint=:zeromean)

Vj = FESpace(
  reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
  conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

Vφ = FESpace(
  reffe=:QLagrangian, order=order-1, valuetype=Float64,
  conformity=:L2, model=model)

U = TrialFESpace(Vu,g_u)
P = TrialFESpace(Vp)
J = TrialFESpace(Vj,g_j)
Φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, J, Φ])

trian = Triangulation(model)
degree = 2*(order)
quad = CellQuadrature(trian,degree)

btrian_j = BoundaryTriangulation(model,"neumann_j")
bquad_j = CellQuadrature(btrian_j,degree)
n = get_normal_vector(btrian_j)

res_Ω(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, f_u)
jac_Ω(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B)
res_Γ(x,y) = InductionlessMHD.dimensionless_conducting_wall(x, y, n, c_w)
jac_Γ(x,dx,y) = InductionlessMHD.dimensionless_conducting_wall(dx, y, n, c_w)

t_Ω = FETerm(res_Ω,jac_Ω,trian,quad)
t_Γ = FETerm(res_Γ,jac_Γ,btrian_j,bquad_j)
op  = FEOperator(X,Y,t_Ω,t_Γ)

nls = NLSolver(;
show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

xh = solve(solver,op)

# if resultsfile != nothing
uh, ph, jh, φh = xh
divj = (∇⋅jh)

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0


# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 20
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
∂p∂z = -Re * K / L^3
u0(x) = analytical_shercliff_u(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

j0(x) = analytical_shercliff_j(side_wall_semilength, hartmann_wall_semilength,
  wall_thickness, wall_conductivity, σ, ν*ρ, ∂p∂z, Ha, number_fourier_sumands, x)

writevtk(trian, "results.vtu",
  cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh, "u"=>u0, "j"=>j0])

# eu_l2, ej_l2 = compute_u_j_errors(uh, jh, u0, j0, trian, quad)
# e_divj = sqrt(sum(integrate(divj*divj,trian,quad)))
#
# @test eu_l2 < 0.3
# @test ej_l2 < 5
# @test e_divj < 1e-12

end #module

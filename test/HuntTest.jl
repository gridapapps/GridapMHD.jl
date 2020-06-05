module HuntTest

using GridapMHD
using GridapMHD.InductionlessMHD

using Gridap
using IncompleteLU
using LineSearches: BackTracking
using Test

# Problem setting
ρ = 1.0
ν = 1.0
σ = 1.0
U0 = 10.0
B0 = 10.0
L = 1.0

Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))
N = Ha^2/Re
K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1))
∂p∂z = -L^3 * K / Re

f_u(x) = VectorValue(0.0,0.0, -∂p∂z) * L/U0^2
g_u = VectorValue(0.0,0.0,0.0)
g_j = VectorValue(0.0,0.0,0.0)
g_φ = 0.0
B = VectorValue(0.0,Ha,0.0)/B0

# Analyical solutions
side_wall_semilength = 1.0
hartmann_wall_semilength = 1.0
number_fourier_sumands = 10
u0(x) = hunt_u(side_wall_semilength, hartmann_wall_semilength, ρ*ν, ∂p∂z,
    Ha, number_fourier_sumands, x)

j0(x) = hunt_j(side_wall_semilength, hartmann_wall_semilength, σ, ρ*ν, ∂p∂z,
    Ha, number_fourier_sumands, x)

# Discretizatoin
order = 2
domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                     sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])*2/sqrt(2)


dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
dirichlet_tags_j = append!(collect(1:20),[25,26])

partition=(5,5,3)
model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false,true))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)

Vu = FESpace(
  reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
  conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

Vp = FESpace(
  reffe=:PLagrangian, order=order-1, valuetype=Float64,
  conformity=:L2, model=model)

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

# Integration
trian = Triangulation(model)
degree = 2*(order)
quad = CellQuadrature(trian,degree)

res(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, f_u)
jac(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B)

t_Ω = FETerm(res,jac,trian,quad)
op  = FEOperator(X,Y,t_Ω)

# Solver
nls = NLSolver(GmresSolver(initially_zero=false,preconditioner=ilu,τ=1e-6);
  show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

xh = solve(solver,op)

# Postprocess
uh, ph, jh, φh = xh
divj = (∇*jh)

# Scale unknowns
uh = uh * U0
jh = jh * σ * B0 * U0

eu = uh - u0
ej = jh - j0

l2(v) = v*v

eu_l2 = sqrt(sum(integrate(l2(eu),trian,quad)))
ej_l2 = sqrt(sum(integrate(l2(ej),trian,quad)))
e_divj = sqrt(sum(integrate(l2(divj),trian,quad)))

@test eu_l2 < 0.0006
@test ej_l2 < 0.008
@test e_divj < 1e-10
end #module

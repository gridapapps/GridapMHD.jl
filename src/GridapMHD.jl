module GridapMHD

using Gridap
using LineSearches: BackTracking

include("PeriodicBC.jl")
using .PeriodicBC

include("Defaults.jl")
using .Defaults
# using Gridap: ∇, divergence

export main

function main(;
              partition::NTuple{3,Int}=(4,4,3),
              map::Function=identity,
              domain::NTuple{6,Float64}=(-1.0,1.0,-1.0,1.0,0.0,0.3),
              periodic_dir::Vector{Int}=[],
              Δt::Float64=1e-4,
              num_time_steps::Int=4,
              maxit::Int=5,
              use_dimensionless_formulation::Bool=true,
              ν::Float64=1.0,
              ρ::Float64=1.0,
              σ::Float64=1.0,
              L::Float64=1.0,
              U::Float64=1.0,
              Re::Float64=1.0,
              Ha::Float64=1.0,
              f::Function=default_f,
              B0::Function=default_B,
              dirichlet_tags::Array{Array{T,1},1}=[collect(1:26),collect(1:26)],
              u0::Function=default_u_ic,
              gx::Function=default_x_dbc,
              ∂x∂n::Function=default_x_nbc,
              write_output::Bool=true,
              output_filename::String="results"
              ) where T

if Δt == 0
  Δt_inv = 0.0
  nt = 1
else
  Δt_inv = 1.0/Δt
  nt = num_time_steps
end
t0 = 0.0
tf = Δt*nt

N = (Ha^2/Re)

order = 2

if length(periodic_dir) > 0
  model = CartesianDiscreteModel(domain,partition,periodic_dir,map)
else
  model = CartesianDiscreteModel(domain,partition,map)
end

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags[1])
add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags[2])

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
    conformity=:L2, model=model, constraint=:zeromean)


trian = Triangulation(model)
degree = 2*(order+1)
quad = CellQuadrature(trian,degree)

gu(x) = gx(x)[1]
gj(x) = gx(x)[2]
U = TrialFESpace(Vu,gu)
P = TrialFESpace(Vp)
j = TrialFESpace(Vj,gj)
φ = TrialFESpace(Vφ)
un = interpolate(Vu, u0)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, j, φ])

neumanntags = setdiff(collect(1:26),dirichlet_tags)
btrian = BoundaryTriangulation(model,neumanntags)
degree = 2*(order-1)
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)


# @law function jxB(x)
#   a = analytical_j(x)
#   b = B(x)
#   return VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])
# end
@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

if use_dimensionless_formulation
  C_ν = (1/Re)
  C_j = N
  C_f = 1/(Re*Re)
  B_0 = 1/Ha
else
  C_ν = ν
  C_j = 1/ρ
  C_f = 1
  B_0 = 1
end

x = get_physical_coordinate(trian)
function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  # uk*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p
  Δt_inv*u*v_u + C_ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p - C_j * vprod(j,B0(x)*B_0)*v_u +
  j*v_j + ∇(φ)*v_j - vprod(u,B0(x)*B_0)*v_j - ∇(v_φ)*j
end

@law conv(u,∇u) = (∇u')*u
@law dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

c(u,v) = inner(v,conv(u,∇(u)))
dc(u,du,v) = inner(v,dconv(du,∇(du),u,∇(u)))

f_u(x) = f(x)[1]
f_p(x) = f(x)[2]
f_j(x) = f(x)[3]
f_φ(x) = f(x)[4]

function l(y)
  v_u, v_p, v_j, v_φ = y
  Δt_inv*un*v_u + v_u*f_u(x)*C_f + v_p*f_p(x) + v_j*f_j(x) + v_φ*f_φ(x)
end

u_nbc(x) = ∂x∂n(x)[1]
p_nbc(x) = ∂x∂n(x)[2]
j_nbc(x) = ∂x∂n(x)[3]
φ_nbc(x) = ∂x∂n(x)[4]

function l_Γ(y)
  v_u, v_p, v_j, v_φ = y
  u_nbc(x) * v_u + p_nbc(x) * v_p + j_nbc(x) * v_j + φ_nbc(x) * v_φ
end

function res(X,Y)
  u   = X[1]
  v_u = Y[1]
  a(X,Y) + c(u, v_u) - l(Y)
end

function jac(X,Y,dX)
  u   = X[1]
  v_u = Y[1]
  du  = dX[1]
  a(dX,Y) + dc(u,du,v_u)
end

t_Ω = FETerm(res,jac,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op  = FEOperator(X,Y,t_Ω)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=maxit)
solver = FESolver(nls)

timeSteps = collect(t0:Δt:tf)

if write_output
  writePVD("results",timeSteps[1:end])
  writevtk(trian, "results/time_"*string(t0)*".vtu",cellfields=["u"=>un])
end

for t in timeSteps[2:end]

  xh = solve(solver,op)
  un, pn, jn, φn = xh

  if (write_output)
    writevtk(trian,"results/time_"*string(t)*".vtu",cellfields=["u"=>un,
                                                                "p"=>pn,
                                                                "j"=>jn,
                                                                "phi"=>φn])
  end

  @show t

  # Update operator
  op = FEOperator(X,Y,t_Ω)
end


end
end # module

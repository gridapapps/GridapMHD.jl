module GridapMHD

using Gridap
using Gridap: ∇, Δ
using LineSearches: BackTracking, Static


include("Defaults.jl")
using .Defaults
# using Gridap: ∇, divergence

export main
export vprod

@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

function main(;
  domain::NTuple{6,Float64}=(-0.5,0.5,-0.5,0.5,-0.5,0.5),
  periodic_dir::Vector=[],
  partition::NTuple{3,Int}=(4,4,3),
  map::Function=identity,
  Δt::Float64=1e-4,
  num_time_steps::Int=4,
  maxit::Int=5,
  order::Int=2,
  dirichlet_tags_u::Vector=collect(1:26),
  dirichlet_tags_j::Vector=collect(1:26),
  use_dimensionless_formulation::Bool=true,
  B::Union{Function,VectorValue}=VectorValue(1.0,1.0,1.0),
  ν::Float64=1.0,
  ρ::Float64=1.0,
  σ::Float64=1.0,
  L::Float64=1.0,
  U::Float64=1.0,
  Re::Float64=1.0,
  Ha::Float64=1.0,
  f_u::Function=default_f_u,
  f_p::Function=default_f_p,
  f_j::Function=default_f_j,
  f_φ::Function=default_f_φ,
  g_u::Function=default_g,
  g_j::Function=default_g,
  ∇u_n::Function=default_∇u_n,
  g_p::Function=default_p,
  g_φ::Function=default_φ,
  u0::Function=default_u_ic,
  write_output::Bool=true,
  output_filename::String="results"
  )

if use_dimensionless_formulation
  C_d = ν
  C_p = 1.0
  C_l = 1/ρ
  C_c = σ
else
  N = Ha^2/Re
  C_d = 1/Re
  C_p = N
  C_l = N
  C_c = 1.0
end


boundary_tags = collect(1:26)
if length(periodic_dir) > 0
  isperiodic = Tuple([(x in periodic_dir) for x in 1:3])
  model = CartesianDiscreteModel(domain,partition;isperiodic=isperiodic,map=map)
  if (1 in periodic_dir) setdiff!(boundary_tags,[25,26]) end
  if (2 in periodic_dir) setdiff!(boundary_tags,[23,24]) end
  if (3 in periodic_dir) setdiff!(boundary_tags,[21,22]) end
else
  model = CartesianDiscreteModel(domain,partition,map)
end

neumann_tags_u = setdiff(boundary_tags,dirichlet_tags_u)
neumann_tags_j = setdiff(boundary_tags,dirichlet_tags_j)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)
add_tag_from_tags!(labels,"neumann_u",neumann_tags_u)
add_tag_from_tags!(labels,"neumann_j",neumann_tags_j)

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet_u")

constraint = :zeromean; if (length(neumann_tags_u) > 0) constraint = nothing end;
Vp = FESpace(
    reffe=:PLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=constraint)

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

constraint = :zeromean; if (length(neumann_tags_j) > 0) constraint = nothing end;
Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=constraint)

U = TrialFESpace(Vu,g_u)
P = TrialFESpace(Vp)
J = TrialFESpace(Vj,g_j)
Φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, J, Φ])

trian = Triangulation(model)
degree = 2*(order)
quad = CellQuadrature(trian,degree)

# uk = interpolate(U,u0)
function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  # (∇(u)'*uk)*v_u +
  C_d*inner(∇(u),∇(v_u)) - C_p*p*(∇*v_u) - C_l*vprod(j,B)*v_u +
  (∇*u)*v_p +
  j*v_j - C_c*φ*(∇*v_j) - C_c*vprod(u,B)*v_j +
  (∇*j)*v_φ
end

function l(Y)
  v_u, v_p, v_j, v_φ = Y

  v_u*f_u + v_p*f_p + v_j*f_j + v_φ*f_φ
end

@law conv(u,∇u) = Re*(∇u')*u
@law dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)
c(x,y) = y[1]*conv(x[1],∇(x[1]))
dc(x,dx,y) = y[1]*dconv(dx[1],∇(dx[1]),x[1],∇(x[1]))

res(x,y) = a(x,y) + c(x,y) - l(y)
jac(x,dx,y) = a(dx,y) + dc(x,dx,y)

t_Ω = FETerm(res,jac,trian,quad)


if (length(neumann_tags_u) > 0)
  btrian_u = BoundaryTriangulation(model,"neumann_u")
  bquad_u = CellQuadrature(btrian_u,degree)
  nb_u = get_normal_vector(btrian_u)

  function l_Γ_u(Y)
    v_u, v_p, v_j, v_φ = Y

    C_d*v_u*(∇u_n) - C_p*(nb_u*v_u)*g_p
  end
  t_Γ_u = FESource(l_Γ_u,btrian_u,bquad_u)
end

if (length(neumann_tags_j) > 0)
  btrian_j = BoundaryTriangulation(model,"neumann_j")
  bquad_j = CellQuadrature(btrian_j,degree)
  nb_j = get_normal_vector(btrian_j)

  function l_Γ_j(Y)
    v_u, v_p, v_j, v_φ = Y

    -C_c*(v_j*nb_j)*g_φ
  end
  t_Γ_j = FESource(l_Γ_j,btrian_j,bquad_j)
end

if (length(neumann_tags_u) > 0) & (length(neumann_tags_j) > 0)
  op  = FEOperator(X,Y,t_Ω, t_Γ_u, t_Γ_j )
elseif (length(neumann_tags_u) > 0)
  op  = FEOperator(X,Y,t_Ω, t_Γ_u )
elseif (length(neumann_tags_j) > 0)
  op  = FEOperator(X,Y,t_Ω, t_Γ_j )
else
  op  = FEOperator(X,Y,t_Ω)
end

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

xh = solve(solver,op)
uh, ph, jh, φh = xh

return xh, trian, quad

end # main
end # module

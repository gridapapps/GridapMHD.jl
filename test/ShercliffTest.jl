module ShercliffTest

include("../src/GridapMHD.jl")
using .GridapMHD

include("../src/Defaults.jl")

using Gridap
using IterativeSolvers
using IncompleteLU
using Preconditioners
using LineSearches: BackTracking
using Test
using Polynomials: fit


ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 10.0
L = 1.0
N = Ha^2/Re
K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1)) # Shercliff's
# K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1)) # Hunt'ss
f_u(x) = VectorValue(0.0,0.0, -L^3 * K / Re)
g_u = VectorValue(0.0,0.0,0.0)
g_j = VectorValue(0.0,0.0,0.0)
B = VectorValue(0.0,1.0,0.0)

# Analyical solutions
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

# Discretizatoin
order = 2
ns = [10]
domain = (-0.5,0.5,-0.5,0.5,0.0,0.1)
map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                     sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])


dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
dirichlet_tags_j = append!(collect(1:20),[23,24,25,26])

##
eu_l2 = Vector{Float64}()
ej_l2 = Vector{Float64}()

for n in ns
  partition=(n,n,3)
  model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false,true),map=map)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
  add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)

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

  U = TrialFESpace(Vu,g_u)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,g_j)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  quad = CellQuadrature(trian,degree)

  function a(X,Y)
    u  , p  , j  , φ   = X
    v_u, v_p, v_j, v_φ = Y

    (1/Re)*inner(∇(u),∇(v_u)) - p*(∇*v_u) - N*vprod(j,B)*v_u +
    (∇*u)*v_p +
    j*v_j - φ*(∇*v_j) - vprod(u,B)*v_j +
    (∇*j)*v_φ
  end

  function l(Y)
    v_u, v_p, v_j, v_φ = Y

    v_u*f_u
  end

  @law conv(u,∇u) = (∇u')*u
  @law dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)
  c(x,y) = y[1]*conv(x[1],∇(x[1]))
  dc(x,dx,y) = y[1]*dconv(dx[1],∇(dx[1]),x[1],∇(x[1]))

  res(x,y) = a(x,y) + c(x,y) - l(y)
  jac(x,dx,y) = a(dx,y) + dc(x,dx,y)

  t_Ω = FETerm(res,jac,trian,quad)
  # t_Ω = AffineFETerm(a,l,trian,quad)

  op  = FEOperator(X,Y,t_Ω)
  # op  = AffineFEOperator(X,Y,t_Ω)

  # A = get_matrix(op)
  # b = get_vector(op)

  # p = ilu(A, τ=1e-6)
  # x = gmres(A b,verbose=true, Pl=p)
  #
  # xh = FEFunction(X,x)

  nls = NLSolver(
    show_trace=true, method=:newton, linesearch=BackTracking())
  solver = FESolver(nls)

  xh = solve(solver,op)

  uh, ph, jh, φh = xh

  eu = uh - u0
  ej = jh - j0
  divj = (∇*jh)

  l2(v) = v*v

  append!(eu_l2, sqrt(sum(integrate(l2(eu),trian,quad))))
  append!(ej_l2, sqrt(sum(integrate(l2(ej),trian,quad))))

  writevtk(trian,"results", nsubcells=order,cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
  "u"=>u0, "j"=>j0, "eu"=>eu, "ej"=>ej, "divj"=>divj])
end

slope_eu_l2 = fit([log(x) for x in ns], [log(y) for y in eu_l2], 1)[end]
slope_ej_l2 = fit([log(x) for x in ns], [log(y) for y in ej_l2], 1)[end]
@show eu_l2
@show ej_l2


end #module

module PeriodicDarcyTests

include("../src/GridapMHD.jl")
using .GridapMHD

using Test
using Gridap
import Gridap: ∇, divergence

function main(u,p,f,g,order_rt,order_p)
  domain = (0,1,0,1,0)
  partition = (20,20)
  order = 0
  model = CartesianDiscreteModel(domain,partition, [1])

  V = FESpace(
  reffe=:RaviartThomas, order=order_rt, valuetype=VectorValue{2,Float64},
  conformity=:Hdiv, model=model,dirichlet_tags=[5,6])
  # conformity=:Hdiv, model=model,dirichlet_tags=append!(collect(1:20),[23,24]))

  Q = FESpace(
  reffe=:QLagrangian, order=order_p, valuetype=Float64,
  conformity=:L2, model=model, constraint=:zeromean)

  U = TrialFESpace(V,u)
  P = TrialFESpace(Q)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  trian = Triangulation(model)
  degree = 2*(order_rt+1)
  quad = CellQuadrature(trian,degree)

  # neumanntags = [26,]
  neumanntags = [7,8]
  btrian = BoundaryTriangulation(model,neumanntags)
  degree = 2*(order_rt+1)
  bquad = CellQuadrature(btrian,degree)
  nb = get_normal_vector(btrian)
  x = get_physical_coordinate(trian)

  function a(x,y)
    u, p = x
    v, q = y
    v*u - p*(∇*v) + q*(∇*u)
  end

  @law f_u(x,v) = v*f(x)
  @law f_p(x,v) = v*g(x)
  function l(y)
    v, q = y
    f_u(x,v) + f_p(x,q)
  end

  function l_Γ(y)
    v, q = y
    -(v*nb)*p
  end

  t_Ω = AffineFETerm(a,l,trian,quad)
  # t_Γ = FESource(l_Γ,btrian,bquad)
  op = AffineFEOperator(X,Y,t_Ω)
  xh = solve(op)
  uh, ph = xh

  writevtk(trian,"results",cellfields=["u"=>uh, "p"=>ph])

  # _u = interpolate(V,u)
  # _p = interpolate(Q,p)
  return trian, quad, uh, ph
end

l2(v) = v*v
tol = 1.0e-9

u(x) = VectorValue(1.0,1.0)
p(x) = 0.0
f(x) = VectorValue(1.0,1.0)
g(x) = 0.0

trian,quad, uh, ph = main(u,p,f,g,0,0)
eu = u - uh
ep = p - ph

eu_l2 = sum(integrate(l2(eu),trian,quad))
ep_l2 = sum(integrate(l2(ep),trian,quad))

@test eu_l2 < tol
@test ep_l2 < tol

##
u(x) = VectorValue(x[1]*(x[1]-1)*(2x[2]-1.0),-x[2]*(x[2]-1.0)*(2x[1]-1.0))
p(x) = x[2]-0.5
f(x) = VectorValue(x[1]*(x[1]-1)*(2x[2]-1.0),-x[2]*(x[2]-1.0)*(2x[1]-1.0)+1.0)
g(x) = 0.0


trian,quad, uh, ph = main(u,p,f,g,1,1)
eu = u - uh
ep = p - ph

eu_l2 = sum(integrate(l2(eu),trian,quad))
ep_l2 = sum(integrate(l2(ep),trian,quad))


@show eu_l2^2
@show ep_l2^2
@test eu_l2 < tol
@test ep_l2 < tol

##

u(x) = VectorValue(2π*cos(2π*x[1])*sin(2π*x[2]), -2π*sin(2π*x[1])*cos(2π*x[2]))
p(x) = sin(2π*x[1])*sin(2π*x[2])


f(x) = VectorValue(4π*cos(2π*x[1])*sin(2π*x[2]),0.0)
g(x) = 0.0

trian,quad, uh, ph = main(u,p,f,g,3,3)
eu = u - uh
ep = p - ph

eu_l2 = sum(integrate(l2(eu),trian,quad))
ep_l2 = sum(integrate(l2(ep),trian,quad))

@test eu_l2 < tol
@test ep_l2 < tol


end # module

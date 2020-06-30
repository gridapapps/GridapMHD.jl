

g_u(x,t) = VectorValue(0.0,0.0,0.0)
g_u(t::Real) = x -> g_u(x,t)
∂tu(t) = x -> VectorValue(0.0,0.0,0.0)
∂tu(x,t) = ∂tu(t)(x)
∂t(::typeof(g_u)) = ∂tu

g_j(x,t) = VectorValue(0.0,0.0,0.0)
g_j(t::Real) = x -> g_j(x,t)
∂tj(t) = x -> VectorValue(0.0,0.0,0.0)
∂tj(x,t) = ∂tj(t)(x)
∂t(::typeof(g_j)) = ∂tu

function transient_duct_flow(;nx::Int=3, ny::Int=3, Re::Float64 = 10.0,
    Ha::Float64 = 10.0, U0::Float64 = Re, B0::Float64 = Ha, L::Float64 = 1.0,
    Δt::Float64 = 0.1, t0::Float64 = 0.0, tF::Float64 = 1.0, θ::Float64 = 1.0,
    resultsfile = nothing)

  N = Ha^2/Re
  K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
  ∂p∂z = -Re * K / L^3

  B = VectorValue(0.0,Ha,0.0)/B0

  f_u(x) = VectorValue(0.0,0.0, -∂p∂z) * L/U0^2

  # Discretization
  order = 2
  domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
  map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                       sign(x[2])*(abs(x[2])*0.5)^0.5, x[3]*sqrt(2)/2)*2/sqrt(2)


  dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
  dirichlet_tags_j = append!(collect(1:20),[23,24,25,26])

  partition=(nx,ny,3)
  model = CartesianDiscreteModel(domain,partition;
    isperiodic=(false,false,true),map=map)

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

  U = TransientTrialFESpace(Vu,g_u)
  P = TransientTrialFESpace(Vp)
  J = TransientTrialFESpace(Vj,g_j)
  Φ = TransientTrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  quad = CellQuadrature(trian,degree)

  res(t,x,xt,y) = xt[1]⋅y[1] +
    InductionlessMHD.dimensionless_residual(x, y, Re, N, B, f_u)
  jac(t,x,xt,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B)
  jac_t(t,x,xt,dxt,y) = dxt[1]⋅y[1]

  t_Ω = FETerm(res,jac,jac_t,trian,quad)
  op  = TransientFEOperator(X,Y,t_Ω)

  uh0 = interpolate(U(0.0),VectorValue(0.0,0.0,0.0))
  ph0 = interpolate(P(0.0),0.0)
  jh0 = interpolate(J(0.0),VectorValue(0.0,0.0,0.0))
  φh0 = interpolate(Φ(0.0),0.0)
  xh0 = Gridap.MultiField.MultiFieldFEFunction(X(0.0),[uh0,ph0,jh0,φh0])

  ls = LUSolver()
  nls = NLSolver(ls;
    show_trace=true, method=:newton, linesearch=BackTracking())
  odes = ThetaMethod(nls,Δt,θ)
  solver = TransientFESolver(odes)

  xh_t = solve(solver,op,xh0,t0,tF)

  if resultsfile != nothing
    startPVD(resultsfile)
    for (xh, t) in xh_t
      uh, ph, jh, φh = xh
      write_timestep(resultsfile,t,trian,
        cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh])
    end
    closePVD(resultsfile)
  end

  (xh_t, trian, quad)
end

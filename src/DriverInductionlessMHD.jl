

function driver_inductionless_MHD(;meshfile=nothing, nx = 4, Re::Float64 = 10.0,
  Ha::Float64 = 10.0, boundary_labels = [], boundary_condition = [],
  boundary_type = [], resultsfile = nothing)

  N = Ha^2/Re
  K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))

  g_u = VectorValue(0.0,0.0,0.0)
  g_j = VectorValue(0.0,0.0,0.0)
  B = VectorValue(0.0,1.0,0.0)

  # Discretization
  order = 2
  if meshfile === nothing
    domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
    map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                        sign(x[2])*(abs(x[2])*0.5)^0.5, x[3]*sqrt(2)/2)*2/sqrt(2)


    dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
    dirichlet_tags_j = append!(collect(1:20),[23,24,25,26])

    partition=(nx,ny,3)
    model = CartesianDiscreteModel(domain,partition;isperiodic=(false,false,true))

    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
    add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)
  else
    model = GmshDiscreteModel(meshfile)

    labels = get_face_labeling(model)
  end


  Vu = FESpace(
      reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
      conformity=:H1, model=model, dirichlet_tags=["side_wall","hartmann_wall"])

  Vp = FESpace(
      reffe=:PLagrangian, order=order-1, valuetype=Float64,
      conformity=:L2, model=model, constraint=:zeromean)

  Vj = FESpace(
      reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
      conformity=:Hdiv, model=model, dirichlet_tags=["side_wall","hartmann_wall"])

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

  res(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, f_u)
  jac(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B)

  t_Ω = FETerm(res,jac,trian,quad)
  op  = FEOperator(X,Y,t_Ω)

  nls = NLSolver(;
    show_trace=true, method=:newton, linesearch=BackTracking())
  solver = FESolver(nls)

  xh = solve(solver,op)

  if resultsfile != nothing
    uh, ph, jh, φh = xh
    writevtk(trian, resultsfile,
      cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh])
  end

  (xh, trian, quad)
end

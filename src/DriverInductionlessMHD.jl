

function driver_inductionless_MHD(;model=nothing, nx = 4, Re::Float64 = 10.0,
  Ha::Float64 = 10.0, fluid_dirichlet_tags = [], fluid_neumann_tags = [],
  magnetic_dirichlet_tags = [], magnetic_neumann_tags = [],
  fluid_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  magnetic_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  fluid_body_force = (x) -> VectorValue(0.0,0.0,0.0),
  constraint_presures::NTuple{2,Bool}=(false,false),
  usegmres = true, precond_tau = 1e-9, linesearch=BackTracking(),
  resultsfile = nothing)

  N = Ha^2/Re
  K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))

  B = VectorValue(0.0,1.0,0.0)

  # Discretization
  order = 2
  if model === nothing
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
  end

  Vu = FESpace(model, ReferenceFE(:Lagrangian,VectorValue{3,Float64},order);
      conformity=:H1, dirichlet_tags=fluid_dirichlet_tags)
  if constraint_presures[1]
    Vp = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:P);
    conformity=:L2, constraint=:zeromean)
  else
    Vp = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:P);
    conformity=:L2)
  end

  Vj = FESpace(model, ReferenceFE(:RaviartThomas,Float64,order-1);
      conformity=:Hdiv, dirichlet_tags=magnetic_dirichlet_tags)
  if constraint_presures[2]
    Vφ = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:Q);
      conformity=:L2, constraint=:zeromean)
  else
    Vφ = FESpace(model, ReferenceFE(:Lagrangian,Float64,order-1,space=:Q);
    conformity=:L2)
  end

  U = TrialFESpace(Vu,fluid_dirichlet_conditions)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,magnetic_dirichlet_conditions)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  dΩ = Measure(trian,degree)

  res(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, fluid_body_force, dΩ)
  jac(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B, dΩ)

  op  = FEOperator(res,jac,X,Y)

  if usegmres
    nls = NLSolver(GmresSolver(preconditioner=ilu,τ=precond_tau);
      show_trace=true, method=:newton, linesearch=linesearch, iterations=10)
  else
    nls = NLSolver(;
      show_trace=true, method=:newton, linesearch=linesearch, iterations=10)
  end
  solver = FESolver(nls)

  xh = solve(solver,op)

  if resultsfile ≠ nothing
    uh, ph, jh, φh = xh
    writevtk(trian, resultsfile,
      cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh])
  end

  (xh, trian, dΩ)
end

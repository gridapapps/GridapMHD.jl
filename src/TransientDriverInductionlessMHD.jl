

function transient_driver_inductionless_MHD(;t0::Float64 = 0.0, tF::Float64 = 1.0,
  model=nothing, nx = 4, Δt::Float64 = 0.1, θ::Float64 = 1.0, Re::Float64 = 10.0,
  Ha = (x) -> 10.0, fluid_dirichlet_tags = [], fluid_neumann_tags = [],
  c_w = 1.0, α=10.0,B = (x) -> VectorValue(0.0,1.0,0.0),
  magnetic_dirichlet_tags = [], magnetic_neumann_tags = [],
  magnetic_non_perfectly_conducting_walls_tag = [],
  fluid_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  magnetic_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  fluid_body_force = (x) -> VectorValue(0.0,0.0,0.0),
  u0 = (x) -> VectorValue(0.0,0.0,0.0), p0 = (x) -> 0.0,
  j0 = (x) -> VectorValue(0.0,0.0,0.0), φ0 = (x) -> 0.0,
  constraint_presures::NTuple{2,Bool}=(false,false), max_nl_it=10,
  usegmres = true, precond_tau = 1e-9, linesearch=BackTracking(),
  resultsfile = nothing, verbosity::Verbosity=Verbosity(1),nsubcells=2)

  if typeof(Ha) == typeof(10.0)
    N = (Ha*Ha)/Re
  else
    N(x) = (Ha(x)*Ha(x))/Re
  end


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
    fluid_dirichlet_tags = ["dirichlet_u"]
    magnetic_dirichlet_tags = ["dirichlet_j"]
  end

  preprocessStepOutput(verbosity,"Preprocess step (0/3): Starting preprocess")

  Vu = FESpace(model, ReferenceFE(lagrangian,VectorValue{3,Float64},order);
      conformity=:H1, dirichlet_tags=fluid_dirichlet_tags)

  if is_n_cube(model.grid_topology.polytopes[1])
    p_conformity = :L2
  else
    p_conformity = :H1
  end

  if constraint_presures[1]
    Vp = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:P);
    conformity=:L2, constraint=:zeromean)
  else
    Vp = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:P);
    conformity=:L2)
  end

  if length(magnetic_dirichlet_tags) == 0
    Vj = FESpace(model, ReferenceFE(raviart_thomas,Float64,order-1);
        conformity=:Hdiv)
  else
    Vj = FESpace(model, ReferenceFE(raviart_thomas,Float64,order-1);
        conformity=:Hdiv, dirichlet_tags=magnetic_dirichlet_tags)
  end

  if constraint_presures[2]
    Vφ = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:Q);
      conformity=:L2, constraint=:zeromean)
  else
    Vφ = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:Q);
    conformity=:L2)
  end

  U = TransientTrialFESpace(Vu,fluid_dirichlet_conditions)
  P = TransientTrialFESpace(Vp)
  J = TransientTrialFESpace(Vj,magnetic_dirichlet_conditions)
  Φ = TransientTrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = TransientMultiFieldFESpace([U, P, J, Φ])

  preprocessStepOutput(verbosity,"Preprocess step (1/3): Built FE spaces")

  trian = Triangulation(model)
  degree = 2*(order)
  dΩ = Measure(trian,degree)

  coord = get_physical_coordinate(trian)
  res_Ω(t,x,xt,y) = ∫(xt[1]⋅y[1])*dΩ +
    InductionlessMHD.dimensionless_residual(x, y, Re, N, B, fluid_body_force,dΩ,coord)
  jac_Ω(t,x,xt,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B,dΩ,coord)
  jac_t(t,x,xt,dxt,y) = ∫(dxt[1]⋅y[1])*dΩ

  if length(magnetic_non_perfectly_conducting_walls_tag) == 0
    op  = TransientFEOperator(res_Ω,jac_Ω,jac_t,X,Y)
  else
    btrian_j = BoundaryTriangulation(model,tags=magnetic_non_perfectly_conducting_walls_tag)
    dΓ = Measure(btrian_j,degree)
    n = get_normal_vector(btrian_j)
    res_Γ(t,x,xt,y) = InductionlessMHD.dimensionless_conducting_wall(x, y, n, c_w, dΓ, α=α)
    jac_Γ(t,x,xt,dx,y) = InductionlessMHD.dimensionless_conducting_wall(dx, y, n, c_w, dΓ, α=α)

    res(t,x,xt,y) = res_Ω(t,x,xt,y) + res_Γ(t,x,xt,y)
    jac(t,x,xt,dx,y) = jac_Ω(t,x,xt,dx,y) + jac_Γ(t,x,xt,dx,y)
    op  = TransientFEOperator(res,jac,jac_t,X,Y)
  end


  uh0 = interpolate(u0,U(0.0))
  ph0 = interpolate(p0,P(0.0))
  jh0 = interpolate(j0,J(0.0))
  φh0 = interpolate(φ0,Φ(0.0))
  xh0 = interpolate([uh0,ph0,jh0,φh0],X(0.0))

  preprocessStepOutput(verbosity,"Preprocess step (2/3): Defined triangulation and formulation")

  if usegmres
    nls = NLSolver(GmresSolver(verbose=linearSolverOutput(verbosity), preconditioner=ilu,τ=precond_tau);
      show_trace=nonlinearSolverOutput(verbosity), method=:newton, linesearch=linesearch, iterations=max_nl_it)
  else
    nls = NLSolver(;
      show_trace=nonlinearSolverOutput(verbosity), method=:newton, linesearch=linesearch, iterations=max_nl_it)
  end
  odes = ThetaMethod(nls,Δt,θ)
  solver = TransientFESolver(odes)

  preprocessStepOutput(verbosity,"Preprocess step (3/3): Configured solver. Starting to solve.")

  xh_t = solve(solver,op,xh0,t0,tF)

  if resultsfile ≠ nothing
    startPVD(resultsfile)
    write_timestep(resultsfile,t0,trian,nsubcells=nsubcells,
        cellfields=["uh"=>uh0, "ph"=>ph0, "jh"=>jh0, "phih"=>φh0])
    it = 0
    timeStepOutput(verbosity,Δt,0.0,it)
    for (xh, t) in xh_t
      it += 1
      timeStepOutput(verbosity,Δt,t,it)

      uh, ph, jh, φh = xh
      write_timestep(resultsfile,t,trian,nsubcells=nsubcells,
        cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "phih"=>φh])
    end
    closePVD(resultsfile)
    return (xh_t, trian, dΩ)

  else
    results = []
    it = 0
    timeStepOutput(verbosity,Δt,0.0,it)
    for (xh, t) in xh_t
      it += 1
      timeStepOutput(verbosity,Δt,t,it)
      results = [results...,xh]
    end
    return (results, trian, dΩ)
  end
end



function driver_inductionless_MHD(;model=nothing, cutgeo=nothing,nx = 4,
  Re::Float64 = 10.0,
  Ha = (x) -> 10.0, fluid_dirichlet_tags = [], fluid_neumann_tags = [],
  c_w = 1.0, α=10.0, B = (x) -> VectorValue(0.0,1.0,0.0),
  magnetic_dirichlet_tags = [], magnetic_neumann_tags = [],
  magnetic_non_perfectly_conducting_walls_tag = [],
  fluid_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  magnetic_dirichlet_conditions = (x) -> VectorValue(0.0,0.0,0.0),
  fluid_body_force = (x) -> VectorValue(0.0,0.0,0.0),
  constraint_presures::NTuple{2,Bool}=(false,false), max_nl_it=10,
  usegmres = true, precond_tau = 1e-9, linesearch=BackTracking(),
  resultsfile = nothing, verbosity::Verbosity=Verbosity(1),nsubcells=1,print_order=2)

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

  grid_topology = get_grid_topology(model)
  if is_n_cube(grid_topology.polytopes[1])
    p_conformity = :L2
  else
    p_conformity = :H1
  end

  if constraint_presures[1]
    Vp = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:P);
    conformity=p_conformity, constraint=:zeromean)
  else
    Vp = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:P);
    conformity=p_conformity)
  end

  if length(magnetic_dirichlet_tags) == 0
    Vj = FESpace(model, ReferenceFE(raviart_thomas,Float64,order-1);
        conformity=:Hdiv)
  else
    Vj = FESpace(model, ReferenceFE(raviart_thomas,Float64,order-1);
        conformity=:Hdiv, dirichlet_tags=magnetic_dirichlet_tags)
  end

  if constraint_presures[2]
    Vφ = FESpace(model, ReferenceFE(lagrangian,Float64,order-1);
      conformity=:L2, constraint=:zeromean)
  else
    Vφ = FESpace(model, ReferenceFE(lagrangian,Float64,order-1);
    conformity=:L2)
  end

  U = TrialFESpace(Vu,fluid_dirichlet_conditions)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,magnetic_dirichlet_conditions)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  preprocessStepOutput(verbosity,"Preprocess step (1/3): Built FE spaces")

  trian = Triangulation(model)
  degree = 2*(order)
  dΩ = Measure(trian,degree)

  coord = get_physical_coordinate(trian)
  res_Ω(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, fluid_body_force, dΩ, coord)
  jac_Ω(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B, dΩ, coord)

  if length(magnetic_non_perfectly_conducting_walls_tag) == 0
    op  = FEOperator(res_Ω,jac_Ω,X,Y)
  else
    if cutgeo === nothing
      btrian_j = BoundaryTriangulation(model,tags=magnetic_non_perfectly_conducting_walls_tag)
    else
      btrian_j = EmbeddedBoundary(cutgeo)
    end
    dΓ = Measure(btrian_j,degree)
    n = get_normal_vector(btrian_j)
    res_Γ(x,y) = InductionlessMHD.dimensionless_conducting_wall(x, y, n, c_w, dΓ, α=α)
    jac_Γ(x,dx,y) = InductionlessMHD.dimensionless_conducting_wall(dx, y, n, c_w, dΓ, α=α)

    res(x,y) = res_Ω(x,y) + res_Γ(x,y)
    jac(x,dx,y) = jac_Ω(x,dx,y) + jac_Γ(x,dx,y)
    op  = FEOperator(res,jac,X,Y)
  end

  preprocessStepOutput(verbosity,"Preprocess step (2/3): Defined triangulation and formulation")

  if usegmres
    nls = NLSolver(GmresSolver(verbose=linearSolverOutput(verbosity), preconditioner=ilu,τ=precond_tau);
      show_trace=nonlinearSolverOutput(verbosity), method=:newton, linesearch=linesearch, iterations=max_nl_it)
  else
    nls = NLSolver(;
      show_trace=nonlinearSolverOutput(verbosity), method=:newton, linesearch=linesearch, iterations=max_nl_it)
  end
  solver = FESolver(nls)

  preprocessStepOutput(verbosity,"Preprocess step (3/3): Configured solver. Starting to solve.")

  xh = solve(solver,op)

  if resultsfile ≠ nothing
    uh, ph, jh, φh = xh
    writevtk(trian, resultsfile, nsubcells=nsubcells, order=print_order,
      cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "phih"=>φh])
  end

  (xh, trian, dΩ)
end



function conductive_thin_wall(;nx::Int=3, ny::Int=3, Re::Float64 = 10.0,
    Ha::Float64 = 10.0, c_w::Float64 = 0.1, α::Float64 = 1.0,
    f_u = nothing, resultsfile = nothing)

  L = 1.0
  N = Ha^2/Re
  K = Ha / (1-0.825*Ha^(-1/2)-Ha^(-1))
  ∂p∂z = -Re * K / L^3

  if f_u === nothing
    _f_u = VectorValue(0.0,0.0, -∂p∂z) * L/Re^2 # Assumed ν=L=1 by default
  else
    _f_u = f_u
  end
  g_u = VectorValue(0.0,0.0,0.0)
  g_j = VectorValue(0.0,0.0,0.0)
  B = VectorValue(0.0,1.0,0.0)

  # Discretization
  order = 2
  domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
  map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,
                       sign(x[2])*(abs(x[2])*0.5)^0.5, x[3]*sqrt(2)/2)*2/sqrt(2)


  dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
  dirichlet_tags_j = append!(collect(1:20),[25,26])
  neumann_j = [23,24]

  partition=(nx,ny,3)
  model = CartesianDiscreteModel(domain,partition;
    isperiodic=(false,false,true), map=map)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
  add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)
  add_tag_from_tags!(labels,"neumann_j",neumann_j)

  Vu = FESpace(model, ReferenceFE(lagrangian,VectorValue{3,Float64},order);
      conformity=:H1, dirichlet_tags="dirichlet_u")

  Vp = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:P);
      conformity=:L2, constraint=:zeromean)

  Vj = FESpace(model, ReferenceFE(raviart_thomas,Float64,order-1);
      conformity=:Hdiv, dirichlet_tags="dirichlet_j")

  Vφ = FESpace(model, ReferenceFE(lagrangian,Float64,order-1,space=:Q);
      conformity=:L2)

  U = TrialFESpace(Vu,g_u)
  P = TrialFESpace(Vp)
  J = TrialFESpace(Vj,g_j)
  Φ = TrialFESpace(Vφ)

  Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
  X = MultiFieldFESpace([U, P, J, Φ])

  trian = Triangulation(model)
  degree = 2*(order)
  dΩ = Measure(trian,degree)

  btrian_j = BoundaryTriangulation(model,tags=["neumann_j"])
  dΓ = Measure(btrian_j,degree)
  n = get_normal_vector(btrian_j)

  res_Ω(x,y) = InductionlessMHD.dimensionless_residual(x, y, Re, N, B, f_u, dΩ)
  jac_Ω(x,dx,y) = InductionlessMHD.dimensionless_jacobian(x, dx, y, Re, N, B, dΩ)
  res_Γ(x,y) = InductionlessMHD.dimensionless_conducting_wall(x, y, n, c_w, dΓ, α=α)
  jac_Γ(x,dx,y) = InductionlessMHD.dimensionless_conducting_wall(dx, y, n, c_w, dΓ, α=α)

  res(x,y) = res_Ω(x,y) + res_Γ(x,y)
  jac(x,dx,y) = jac_Ω(x,dx,y) + jac_Γ(x,dx,y)
  op  = FEOperator(res,jac,X,Y)

  nls = NLSolver(;
    show_trace=true, method=:newton, linesearch=BackTracking())
  solver = FESolver(nls)

  xh = solve(solver,op)

  if resultsfile != nothing
    uh, ph, jh, φh = xh
    writevtk(trian, resultsfile,
      cellfields=["uh"=>uh, "ph"=>ph, "jh"=>jh, "φh"=>φh])
  end

  (xh, trian, dΩ)
end

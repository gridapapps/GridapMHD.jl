
function Li2019Solver(op::FEOperator,params)
  U, V = get_trial(op), get_test(op)
  U_u, U_p, U_j, U_φ = U
  V_u, V_p, V_j, V_φ = V

  if isa(params[:ζ],Nothing)
    params[:ζ] = 0.0
  end

  # Preconditioner
  k  = params[:fespaces][:k]
  γ  = params[:fluid][:γ]
  Ωf = _interior(params[:model],params[:fluid][:domain])
  dΩ = Measure(Ωf,2*k)
  Dj = _Dj(U_j,V_j,Ωf,dΩ,params)
  Fk = _Fk(zero(U_u),U_u,V_u,Ωf,dΩ,params)
  Δp = _p_laplacian(U_p,V_p,Ωf,dΩ,params)
  Ip = assemble_matrix((p,v_p) -> ∫(p*v_p)*dΩ,V_p,V_p)
  Iφ = assemble_matrix((φ,v_φ) -> ∫(-γ*φ*v_φ)*dΩ ,U_φ,V_φ)

  block_solvers = map(s -> get_block_solver(Val(s)),params[:solver][:block_solvers])
  block_mats    = [Dj,Fk,Δp,Ip,Iφ]
  P = Li2019_Preconditioner(block_solvers...,block_mats...,params)
  #test_preconditioner(op,P)

  # Linear Solver
  m = params[:solver][:niter]
  l_solver = GMRESSolver(m,P;tol=1e-8,verbose=i_am_main(get_parts(U)))

  # Nonlinear Solver
  nlsolver = NewtonRaphsonSolver(l_solver,1e-5,10)
  return nlsolver
end

function _Dj(U_j,V_j,Ω,dΩ,params)
  fluid = params[:fluid]
  γ  = fluid[:γ]
  k = params[:fespaces][:k]

  params_thin_wall = []
  bcs = params[:bcs]
  for i in 1:length(bcs[:thin_wall])
    τ   = bcs[:thin_wall][i][:τ]
    cw  = bcs[:thin_wall][i][:cw]
    jw  = bcs[:thin_wall][i][:jw]
    Γ   = _boundary(params[:model],bcs[:thin_wall][i][:domain])
    dΓ  = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_thin_wall,(τ,cw,jw,n_Γ,dΓ))
  end

  function a_j(j,v_j) 
    r = ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ 
    for p in params_thin_wall
      τ,cw,jw,n_Γ,dΓ = p
      r += ∫(τ*(v_j⋅n_Γ)⋅(j⋅n_Γ) + cw*(v_j⋅n_Γ)⋅(n_Γ⋅(∇(j)⋅n_Γ)))*dΓ
    end
    return r
  end
  return assemble_matrix(a_j,U_j,V_j)
end

function _Fk(u,U_u,V_u,Ω,dΩ,params)
  fluid = params[:fluid]
  α  = fluid[:α]
  β  = fluid[:β]
  γ  = fluid[:γ]
  B  = fluid[:B]
  
  conv(u,∇u)  = (∇u')⋅u
  a_fk(u,du,dv) = ∫(β*(∇(du)⊙∇(dv)) + α*dv⋅((conv∘(u,∇(du))) + (conv∘(du,∇(u)))) + γ⋅(du×B)⋅(dv×B)) * dΩ
  return assemble_matrix((du,dv) -> a_fk(u,du,dv),U_u,V_u)
end

function _p_laplacian(U_p,V_p,Ω,dΩ,params)
  p_conformity = params[:fespaces][:p_conformity]
  if p_conformity == :H1
    return assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩ ,U_p,V_p)
  elseif p_conformity == :L2
    # TODO: This is only ok if we do not have a solid/fuid split of the mesh...
    model = params[:model]
    Γ  = Boundary(model)
    Λ  = Skeleton(model)

    k  = params[:fespaces][:k]
    dΓ = Measure(Γ,2*k)
    dΛ = Measure(Λ,2*k)

    n_Γ = get_normal_vector(Γ)
    n_Λ = get_normal_vector(Λ)

    h_e_Λ = get_edge_measures(Λ,dΛ)
    h_e_Γ = get_edge_measures(Γ,dΓ)

    β = 100.0
    aΛ(u,v) = ∫(-jump(u⋅n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ))*dΛ + ∫(β/h_e_Λ*jump(u⋅n_Λ)⋅jump(v⋅n_Λ))*dΛ
    aΓ(u,v) = ∫(-(∇(u)⋅n_Γ)⋅v - u⋅(∇(v)⋅n_Γ))*dΓ + ∫(β/h_e_Γ*(u⋅n_Γ)⋅(v⋅n_Γ))*dΓ

    ap(p,v_p) = ∫(∇(p)⋅∇(v_p))*dΩ + aΛ(p,v_p) + aΓ(p,v_p)
    return assemble_matrix(ap,U_p,V_p)
  else
    error("Unknown p_conformity: $p_conformity")
  end
end

get_edge_measures(Ω::Triangulation,dΩ) = CellField(get_array(∫(1)dΩ),Ω)
function get_edge_measures(Ω::GridapDistributed.DistributedTriangulation,dΩ) 
  return CellField(map(get_array,local_views(∫(1)*dΩ)),Ω)
end

function Gridap.Algebra.solve!(x::AbstractVector,nls::NewtonRaphsonSolver,op::NonlinearOperator,cache::Nothing)
  b  = residual(op, x)
  A  = jacobian(op, x)
  dx = allocate_col_vector(A)
  ns = numerical_setup(symbolic_setup(nls.ls, A), A)

  Gridap.Algebra._solve_nr!(x,A,b,dx,ns,nls,op)
  return Gridap.Algebra.NewtonRaphsonCache(A,b,dx,ns)
end

# TODO: This is copied from main... should be in a common place
function _interior(model,domain::Union{Gridap.DiscreteModel,GridapDistributed.DistributedDiscreteModel})
  return Interior(domain)
end
function _interior(model,domain::Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation})
  return domain
end
function _interior(model,domain)
  return Interior(model,tags=domain)
end

function _boundary(model,domain::Union{Gridap.Triangulation,GridapDistributed.DistributedTriangulation})
  return domain
end
function _boundary(model,domain)
  Boundary(model,tags=domain)
end

function test_preconditioner(op,P)
  parts = GridapDistributed.get_parts(get_trial(op))
  i_am_main(parts) && println(repeat("*",50))
  i_am_main(parts) && println("Testing Preconditioner: ")
  xh = zero(get_trial(op))
  A = jacobian(op,xh)
  b = residual(op,xh)

  # Test blocks
  block_solvers = [P.Dj_solver,P.Fk_solver,P.Δp_solver,P.Ip_solver,P.Iφ_solver]
  block_mats    = [P.Dj,P.Fk,P.Δp,P.Ip,P.Iφ]
  test_block_solvers(parts,block_solvers,block_mats)

  # Test preconditioner
  i_am_main(parts) && println(repeat("=",50))
  ns = numerical_setup(symbolic_setup(P,A),A)
  x = allocate_col_vector(A)
  solve!(x,ns,b)
  e = norm(x)
  i_am_main(parts) && println("P error: ",e)
  GridapPETSc.gridap_petsc_gc()
end

function test_block_solvers(parts,block_solvers,block_mats)
  Dj_s,Fk_s,Δp_s,Ip_s,Iφ_s = block_solvers
  Dj,Fk,Δp,Ip,Iφ = block_mats

  function test_solver(s,m)
    ns = numerical_setup(symbolic_setup(s,m),m)
    x = allocate_col_vector(m)
    b = allocate_col_vector(m)
    fill!(b,1.0)
    x = solve!(x,ns,b)
    return norm(b - m*x)
  end

  # Test Dj solver
  i_am_main(parts) && println(repeat("=",50))
  e = test_solver(Dj_s,Dj)
  i_am_main(parts) && println(" > Dj error:",e)

  # Test Fk solver
  i_am_main(parts) && println(repeat("=",50))
  e = test_solver(Fk_s,Fk)
  i_am_main(parts) && println(" > Fk error:",e)

  # Test Δp solver
  i_am_main(parts) && println(repeat("=",50))
  e = test_solver(Δp_s,Δp)
  i_am_main(parts) && println(" > Δp error:",e)

  # Test Ip solver
  i_am_main(parts) && println(repeat("=",50))
  e = test_solver(Ip_s,Ip)
  i_am_main(parts) && println(" > Ip error:",e)

  # Test Iφ solver
  i_am_main(parts) && println(repeat("=",50))
  e = test_solver(Iφ_s,Iφ)
  i_am_main(parts) && println(" > Iφ error:",e)
end

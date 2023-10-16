
function Li2019Solver(op::FEOperator,params)
  U, V = get_trial(op), get_test(op)
  U_u, U_p, U_j, U_φ = U
  V_u, V_p, V_j, V_φ = V

  nl_rtol = params[:solver][:rtol]
  l_rtol  = nl_rtol/100.0

  if isa(params[:ζ],Nothing)
    params[:ζ] = 0.0
  end

  # Preconditioner
  k  = params[:fespaces][:k]
  γ  = params[:fluid][:γ]
  Ωf = _interior(params[:model],params[:fluid][:domain])
  dΩ = Measure(Ωf,2*k)
  a_Dj = _Dj(dΩ,params)
  a_Fk = _Fk(dΩ,params)
  #a_Δp = _p_laplacian(dΩ,params)
  a_Ip(p,v_p) = ∫(p*v_p)*dΩ
  a_Iφ(φ,v_φ) = ∫(-φ*v_φ)*dΩ

  block_solvers   = map(s -> get_block_solver(Val(s)),params[:solver][:block_solvers])
  block_weakforms = [a_Dj,a_Fk,a_Ip,a_Iφ]
  P = Li2019_Preconditioner(op,block_solvers,block_weakforms,params)
  #test_preconditioner(op,P)

  # Linear Solver
  m = params[:solver][:niter]
  l_solver = FGMRESSolver(m,P;rtol=l_rtol,atol=1e-14,verbose=i_am_main(get_parts(U)))

  # Nonlinear Solver
  nlsolver = NewtonRaphsonSolver(l_solver,nl_rtol,10)
  return nlsolver
end

function _Dj(dΩ,params)
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
    r = ∫(j⋅v_j + (∇⋅j)⋅(∇⋅v_j))*dΩ 
    for p in params_thin_wall
      τ,cw,jw,n_Γ,dΓ = p
      r += ∫(τ*(v_j⋅n_Γ)⋅(j⋅n_Γ) + cw*(v_j⋅n_Γ)⋅(n_Γ⋅(∇(j)⋅n_Γ)))*dΓ
    end
    return r
  end
  return a_j
end

function _Fk(dΩ,params)
  fluid = params[:fluid]
  α  = fluid[:α]
  β  = fluid[:β]
  γ  = fluid[:γ]
  B  = fluid[:B]
  
  conv(u,∇u)  = (∇u')⋅u
  function a_fk(u,du,dv) 
    r = ∫(β*(∇(du)⊙∇(dv)) + α*dv⋅((conv∘(u,∇(du))) + (conv∘(du,∇(u)))) + γ⋅(du×B)⋅(dv×B)) * dΩ
    return r
  end
  return a_fk
end

function _p_laplacian(dΩ,params)
  p_conformity = params[:fespaces][:p_conformity]

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
  function a_Δp(u,v)
    r = ∫(∇(u)⋅∇(v))*dΩ
    if p_conformity == :L2
      r += ∫(-jump(u⋅n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ) + β/h_e_Λ*jump(u⋅n_Λ)⋅jump(v⋅n_Λ))*dΛ
      r += ∫(-(∇(u)⋅n_Γ)⋅v - u⋅(∇(v)⋅n_Γ) + β/h_e_Γ*(u⋅n_Γ)⋅(v⋅n_Γ))*dΓ
    end
    return r
  end
  return a_Δp
end

get_edge_measures(Ω::Triangulation,dΩ) = sqrt∘CellField(get_array(∫(1)dΩ),Ω)
function get_edge_measures(Ω::GridapDistributed.DistributedTriangulation,dΩ)
  return sqrt∘CellField(map(get_array,local_views(∫(1)*dΩ)),Ω)
end

function Gridap.Algebra.numerical_setup!(
    ns::GridapSolvers.LinearSolvers.FGMRESNumericalSetup,
    A::AbstractMatrix,
    x::AbstractVector)
  numerical_setup!(ns.Pl_ns,A,x)
  ns.A = A
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

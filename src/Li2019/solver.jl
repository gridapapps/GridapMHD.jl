
function Li2019Solver(op,params)
  U, V = get_trial(op), get_test(op)
  U_u, U_p, U_j, U_φ = U
  V_u, V_p, V_j, V_φ = V

  # Preconditioner
  k  = params[:k]
  γ  = params[:fluid][:γ]
  Ωf = _interior(params[:model],params[:fluid][:domain])
  dΩ = Measure(Ωf,2*k)
  Dj = assemble_matrix((j,v_j) -> ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ ,U_j,V_j)
  Ij = assemble_matrix((j,v_j) -> ∫(j⋅v_j)*dΩ ,U_j,V_j)
  Δp = assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩ ,U_p,V_p)
  Ip = assemble_matrix((p,v_p) -> ∫(p*v_p)*dΩ,V_p,V_p)
  Iφ = assemble_matrix((φ,v_φ) -> ∫(-γ*φ*v_φ)*dΩ ,U_φ,V_φ)

  block_solvers = map(s -> get_block_solver(Val(s)),params[:solver][:block_solvers])
  block_mats    = [Dj,Δp,Ip,Ij,Iφ]
  P = Li2019_Preconditioner(block_solvers...,block_mats...,params)

  # Linear Solver
  l_solver = GMRESSolver(150,P,1e-8)

  # Nonlinear Solver
  nlsolver = NewtonRaphsonSolver(l_solver,1e-5,10)
  return nlsolver
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
  e = test_solver(Dj_s,Dj)
  i_am_main(parts) && println("Dj error:",e)

  # Test Fk solver
  e = test_solver(Fk_s,Fk)
  i_am_main(parts) && println("Fk error:",e)

  # Test Δp solver
  e = test_solver(Δp_s,Δp)
  i_am_main(parts) && println("Δp error:",e)

  # Test Ip solver
  e = test_solver(Ip_s,Ip)
  i_am_main(parts) && println("Ip error:",e)

  # Test Iφ solver
  e = test_solver(Iφ_s,Iφ)
  i_am_main(parts) && println("Iφ error:",e)
end



struct Li2019_Preconditioner <: Gridap.Algebra.LinearSolver
  op
  block_solvers
  block_weakforms
  params
end

struct LI2019_SS <: Gridap.Algebra.SymbolicSetup
  solver::Li2019_Preconditioner
end

function Gridap.Algebra.symbolic_setup(solver::Li2019_Preconditioner, A::AbstractBlockMatrix)
  return LI2019_SS(solver)
end

mutable struct LI2019_NS <: Gridap.Algebra.NumericalSetup
  solver
  Dj_ns
  Fk_ns
  Δp_ns
  Ip_ns
  Iφ_ns
  sysmat
  caches
end

function allocate_caches(solver::Li2019_Preconditioner,A::AbstractBlockMatrix)
  du = allocate_col_vector(A[Block(1,1)])
  dp = allocate_col_vector(A[Block(2,2)])
  dj = allocate_col_vector(A[Block(3,3)])
  dφ = allocate_col_vector(A[Block(4,4)])
  return du, dp, dj, dφ
end

function Gridap.Algebra.numerical_setup(ss::LI2019_SS, A::AbstractBlockMatrix)
  solver = ss.solver

  U = get_trial(solver.op); U_u, U_p, U_j, U_φ = U
  V = get_test(solver.op);  V_u, V_p, V_j, V_φ = V
  a_Dj,a_Fk,a_Δp,a_Ip,a_Iφ = solver.block_weakforms
  Dj_solver,Fk_solver,Δp_solver,Ip_solver,Iφ_solver = solver.block_solvers

  u  = zero(U_u)
  Dj = assemble_matrix(a_Dj,U_j,V_j)
  Fk = assemble_matrix((du,dv) -> a_Fk(u,du,dv),U_u,V_u)
  Δp = assemble_matrix(a_Δp,U_p,V_p)
  Ip = assemble_matrix(a_Ip,U_p,V_p)
  Iφ = assemble_matrix(a_Iφ,U_φ,V_φ)

  Dj_ns = numerical_setup(symbolic_setup(Dj_solver,Dj),Dj)
  Fk_ns = numerical_setup(symbolic_setup(Fk_solver,Fk),Fk)
  Δp_ns = numerical_setup(symbolic_setup(Δp_solver,Δp),Δp)
  Ip_ns = numerical_setup(symbolic_setup(Ip_solver,Ip),Ip)
  Iφ_ns = numerical_setup(symbolic_setup(Iφ_solver,Iφ),Iφ)
  cache = allocate_caches(solver,A)
  return LI2019_NS(ss.solver,Dj_ns,Fk_ns,Δp_ns,Ip_ns,Iφ_ns,A,cache)
end

function Gridap.Algebra.numerical_setup!(ns::LI2019_NS, A::AbstractBlockMatrix)
  Fk = A[Block(1,1)]
  numerical_setup!(ns.Fk_ns,Fk)
  ns.sysmat = A
  return ns
end

function Gridap.Algebra.numerical_setup!(ns::LI2019_NS, A::AbstractBlockMatrix, x::AbstractBlockVector)
  solver = ns.solver
  U = get_trial(solver.op); U_u, _, _, _ = U
  V = get_test(solver.op);  V_u, _, _, _ = V
  _ ,a_Fk, _, _, _ = solver.block_weakforms

  u  = FEFunction(U_u,x[Block(1)])
  Fk = assemble_matrix((du,dv) -> a_Fk(u,du,dv),U_u,V_u)
  ns.Fk_ns  = numerical_setup!(ns.Fk_ns,Fk)
  ns.sysmat = A
  return ns
end

# Follows Algorithm 4.1 in (Li,2019)
function Gridap.Algebra.solve!(x::AbstractBlockVector,ns::LI2019_NS,b::AbstractBlockVector)
  sysmat, caches, params = ns.sysmat, ns.caches, ns.solver.params
  fluid = params[:fluid]; ζ = params[:ζ]; iRe = fluid[:β]
  κ = fluid[:γ]; α1 = ζ + iRe;

  bu, bp, bj, bφ = blocks(b)
  u , p , j , φ  = blocks(x)
  du, dp, dj, dφ = caches

  # Solve for p
  solve!(dp,ns.Δp_ns,bp)
  solve!(p,ns.Ip_ns,bp)
  p .= -α1 .* p .- dp

  #  Solve for φ
  #dφ .= -bφ
  solve!(φ,ns.Iφ_ns,bφ)

  # Solve for u
  copy!(du,bu)
  mul!(du,sysmat[Block(1,2)],p,-1.0,1.0) # du = bu - Aup * p
  solve!(u,ns.Fk_ns,du)                  # u = Fu \ (bu - Aup * p)

  # Solve for j
  copy!(dj,bj)
  mul!(dj,sysmat[Block(3,1)],u,-2.0,1.0) # dj = bj - 2.0 * Aju * u
  mul!(dj,sysmat[Block(3,4)],φ,-2.0,1.0) # dj = bj - 2.0 * Aju * u - 2.0 * Ajφ * φ
  solve!(j,ns.Dj_ns,dj)                  # j = Dj \ (bj - 2.0 * Aju * u - 2.0 * Ajφ * φ)
  return x
end
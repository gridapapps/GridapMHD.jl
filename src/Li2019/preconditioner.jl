
struct Li2019_Preconditioner <: Gridap.Algebra.LinearSolver
  Dj_solver
  Fk_solver
  Δp_solver
  Ip_solver
  Iφ_solver
  Dj
  Fk
  Δp
  Ip
  Iφ
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

  #Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  #Fk = Fu# - (1.0/κ^2) * Kᵗ * K

  Dj_ns = numerical_setup(symbolic_setup(solver.Dj_solver,solver.Dj),solver.Dj)
  Fk_ns = numerical_setup(symbolic_setup(solver.Fk_solver,solver.Fk),solver.Fk)
  Δp_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Δp),solver.Δp)
  Ip_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Ip),solver.Ip)
  Iφ_ns = numerical_setup(symbolic_setup(solver.Iφ_solver,solver.Iφ),solver.Iφ)
  cache = allocate_caches(solver,A)
  return LI2019_NS(ss.solver,Dj_ns,Fk_ns,Δp_ns,Ip_ns,Iφ_ns,A,cache)
end

function Gridap.Algebra.numerical_setup!(ns::LI2019_NS, A::AbstractBlockMatrix)
  solver = ns.solver

  #! Pattern of matrix changes, so we need to recompute everything.
  # This will get fixed when we are using iterative solvers for Fk
  Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  Fk = Fu# - (1.0/κ^2) * Kᵗ * K
  #numerical_setup!(ns.Fk_ns,Fk)

  ns.Fk_ns  = numerical_setup(symbolic_setup(solver.Fk_solver,Fk),Fk)
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
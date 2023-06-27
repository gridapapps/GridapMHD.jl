
using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra

using LinearAlgebra, BlockArrays

using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector

using GridapMHD: _hunt, add_default_params, _fluid_mesh, weak_form, _find_backend, p_conformity, _interior

function Gridap.Algebra._check_convergence(nls,b,m0)
  m = Gridap.Algebra._inf_norm(b)
  println(">>>>>>>>>>>>>>>>>>>> Nonlinear Abs Error = $m")
  m < nls.tol * m0
end

"""
  This preconditioner is based on [(Li,2019)](https://doi.org/10.1137/19M1260372)
"""
struct MHDBlockPreconditioner <: Gridap.Algebra.LinearSolver
  Dj_solver
  Fk_solver
  Δp_solver
  Ip_solver
  Iφ_solver
  Dj
  Δp
  Ip
  Ij
  Iφ
  params
end

struct MHDBlockPreconditionerSS <: Gridap.Algebra.SymbolicSetup
  solver
end

function Gridap.Algebra.symbolic_setup(solver::MHDBlockPreconditioner, A::AbstractMatrix)
  return MHDBlockPreconditionerSS(solver)
end

mutable struct MHDBlockPreconditionerNS <: Gridap.Algebra.NumericalSetup
  solver
  Dj_ns
  Fk_ns
  Δp_ns
  Ip_ns
  Iφ_ns
  sysmat
  caches
end

function allocate_preconditioner_caches(solver,A::BlockMatrix)
  du = allocate_col_vector(A[Block(1,1)])
  dp = allocate_col_vector(A[Block(2,2)])
  dj = allocate_col_vector(A[Block(3,3)])
  dφ = allocate_col_vector(A[Block(4,4)])
  return du, dp, dj, dφ
end

function Gridap.Algebra.numerical_setup(ss::MHDBlockPreconditionerSS, A::BlockMatrix)
  solver = ss.solver

  Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  Fk = Fu - (1.0/κ) * Kᵗ * K

  Dj_ns = numerical_setup(symbolic_setup(solver.Dj_solver,solver.Dj),solver.Dj)
  Fk_ns = numerical_setup(symbolic_setup(solver.Fk_solver,Fk),Fk)
  Δp_ns = nothing#numerical_setup(symbolic_setup(solver.Δp_solver,solver.Δp),solver.Δp)
  Ip_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Ip),solver.Ip)
  Iφ_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Iφ),solver.Iφ)
  caches = allocate_preconditioner_caches(solver,A)
  return MHDBlockPreconditionerNS(ss.solver,Dj_ns,Fk_ns,Δp_ns,Ip_ns,Iφ_ns,A,caches)
end

function Gridap.Algebra.numerical_setup!(ns::MHDBlockPreconditionerNS, A::BlockMatrix)
  solver = ns.solver

  #! Pattern of matrix changes, so we need to recompute everything.
  # This will get fixed when we are using iterative solvers for Fk
  Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  Fk = Fu - (1.0/κ) * Kᵗ * K
  # numerical_setup!(ns.Fk_ns,Fk)

  ns.Fk_ns  = numerical_setup(symbolic_setup(solver.Fk_solver,Fk),Fk)
  ns.sysmat = A

  return ns
end

# Follows Algorithm 4.1 in (Li,2019)
function Gridap.Algebra.solve!(x::BlockVector,ns::MHDBlockPreconditionerNS,b::BlockVector)
  sysmat, caches, params = ns.sysmat, ns.caches, ns.solver.params
  fluid = params[:fluid]; ζ = params[:ζ]; iRe = fluid[:β]
  κ = fluid[:γ]; α1 = ζ + iRe;

  bu, bp, bj, bφ = blocks(b)
  u, p, j, φ = blocks(x)
  du, dp, dj, dφ = caches

  # Solve for p
  #solve!(p,ns.Δp_ns,bp)
  solve!(dp,ns.Ip_ns,bp)
  p .= -α1 .* dp #.- p #! Time-dependent case:  .+ (2.0/τ) .* p; See Lemma A.4 in paper.

  #  Solve for φ
  dφ .= -bφ
  solve!(φ,ns.Iφ_ns,dφ)

  # Solve for u
  copy!(du,bu); mul!(du,sysmat[Block(1,2)],p,-1.0,1.0) # du = bu - Aup * p
  solve!(u,ns.Fk_ns,du) # u = Fu \ (bu - Aup * p)

  # Solve for j
  copy!(dj,bj)
  mul!(dj,sysmat[Block(3,1)],u,-2.0,1.0) # dj = bj - 2.0 * Aju * u
  mul!(dj,sysmat[Block(3,4)],φ,-2.0,1.0) # dj = bj - 2.0 * Aju * u - 2.0 * Ajφ * φ
  solve!(j,ns.Dj_ns,dj) # j = Dj \ (bj - 2.0 * Aju * u - 2.0 * Ajφ * φ)

  return x
end

function hunt(;
  backend=nothing,
  np=nothing,
  parts=nothing,
  title = "hunt",
  path=".",
  kwargs...)

  @assert parts === nothing
  if backend === nothing
    @assert np === nothing
    return _hunt(;title=title,path=path,kwargs...)
  else
    @assert backend !== nothing
    return with_backend(_find_backend(backend),(np...,1)) do _parts
      _hunt(;parts=_parts,title=_title,path=path,kwargs...)
    end
  end
end

_params = hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=false,
  solver=:block_gmres,
)

params = add_default_params(_params)

# ReferenceFEs
k = params[:k]
T = Float64
model = params[:model]
D = num_cell_dims(model)
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)
reffe_φ = ReferenceFE(lagrangian,T,k-1)

# Test spaces
Ωf  = _fluid_mesh(model,params[:fluid][:domain])
V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
V_p = TestFESpace(Ωf,reffe_p;conformity=p_conformity(Ωf))
V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V   = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=BlockMultiFieldStyle())

# Trial spaces

z = zero(VectorValue{D,Float64})
u_bc = params[:bcs][:u][:values]
j_bc = params[:bcs][:j][:values]
U_u  = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
U_j  = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
U_p  = TrialFESpace(V_p)
U_φ  = TrialFESpace(V_φ)
U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=BlockMultiFieldStyle())

# Weak form

params[:ζ] = 100.0
res, jac = weak_form(params,k)
Tm = params[:matrix_type]
Tv = params[:vector_type]
assem = SparseMatrixAssembler(Tm,Tv,U,V)
op    = FEOperator(res,jac,U,V,assem)
al_op = Gridap.FESpaces.get_algebraic_operator(op)

# Preconditioner

Ω = Triangulation(params[:model])
dΩ = Measure(Ω,2*k)

Dj = assemble_matrix((j,v_j) -> ∫(j⋅v_j + (∇⋅j)⋅(∇⋅v_j))*dΩ ,U_j,V_j)
Ij = assemble_matrix((j,v_j) -> ∫(j⋅v_j)*dΩ ,U_j,V_j)
Δp = assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩ ,U_p,V_p)
Ip = assemble_matrix((p,v_p) -> ∫(p*v_p)*dΩ,V_p,V_p)
Iφ = assemble_matrix((φ,v_φ) -> ∫(φ*v_φ)*dΩ ,U_φ,V_φ)

Dj_solver = LUSolver()
Fk_solver = LUSolver()
Δp_solver = LUSolver() # Not used for now since Δp is singular
Ip_solver = LUSolver()
Iφ_solver = LUSolver()

block_solvers = [Dj_solver,Fk_solver,Δp_solver,Ip_solver,Iφ_solver]
block_mats = [Dj,Δp,Ip,Ij,Iφ]
P = MHDBlockPreconditioner(block_solvers...,block_mats...,params)

sysmat_solver = GMRESSolver(300,P,1e-8)
sysmat_ns = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)


# Gridap's Newton-Raphson solver
xh = zero(U)
sysvec = residual(op,xh)
sysmat = jacobian(op,xh)

x  = allocate_col_vector(sysmat)
dx = allocate_col_vector(sysmat)

nlsolver = NewtonRaphsonSolver(sysmat_solver,1e-6,100)
nlsolver_cache = Gridap.Algebra.NewtonRaphsonCache(sysmat,sysvec,dx,sysmat_ns)
solve!(x,nlsolver,al_op,nlsolver_cache)

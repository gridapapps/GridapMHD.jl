
using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra

using LinearAlgebra, BlockArrays

using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector

using GridapMHD: _hunt, add_default_params, _fluid_mesh, weak_form, _find_backend, p_conformity

struct MHDBlockPreconditioner <: Gridap.Algebra.LinearSolver
  Dj_solver
  Fu_solver
  Δp_solver
  Dj
  Δp
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
  Fu_ns
  Δp_ns
  sysmat
  caches
end

function allocate_preconditioner_caches(solver,A::BlockMatrix)
  du = allocate_col_vector(A[Block(1,1)])
  dj = allocate_col_vector(A[Block(3,3)])
  return du, dj
end

function Gridap.Algebra.numerical_setup(ss::MHDBlockPreconditionerSS, A::BlockMatrix)
  solver = ss.solver
  Fu = A[Block(1,1)]
  Dj_ns = numerical_setup(symbolic_setup(solver.Dj_solver,solver.Dj),solver.Dj)
  Fu_ns = numerical_setup(symbolic_setup(solver.Fu_solver,Fu),Fu)
  Δp_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Δp),solver.Δp)
  caches = allocate_preconditioner_caches(solver,A)
  return MHDBlockPreconditionerNS(ss.solver,Dj_ns,Fu_ns,Δp_ns,A,caches)
end

function Gridap.Algebra.numerical_setup!(ns::MHDBlockPreconditionerNS, A::BlockMatrix)
  numerical_setup!(ns.Fu_ns,A[Block(1,1)])
  ns.sysmat = A
end

function Gridap.Algebra.solve!(x::BlockVector,ns::MHDBlockPreconditionerNS,b::BlockVector)
  sysmat, caches, params = ns.sysmat, ns.caches, ns.solver.params
  fluid = params[:fluid]; α = fluid[:α]; iRe = fluid[:β]
  κ = 1.0; α1 = α + iRe; τ = 1.0

  bu, bp, bj, bφ = blocks(b)
  u, p, j, φ = blocks(x)
  du, dj = caches

  # Solve for φ and p
  φ .= (-1.0/κ) .* bφ
  solve!(p,ns.Δp_ns,bp); p .= -α1 .* bp + (2.0/τ) .* p;

  # Solve for u
  copy!(du,bu); mul!(du,sysmat[Block(1,2)],p,-1.0,1.0) # du = bu - Aup * p
  solve!(u,ns.Fu_ns,du) # u = Fu \ (bu - Aup * p)

  # Solve for j
  copy!(dj,bj)
  mul!(dj,sysmat[Block(3,1)],u,-2.0,1.0) # dj = bj - 2.0 * Aju * u
  mul!(dj,sysmat[Block(3,4)],φ,-2.0,1.0) # dj = bj - 2.0 * Aju * u - 2.0 * Ajφ * φ
  solve!(j,ns.Dj_ns,dj) # j = Dj \ (bj - 2.0 * Aju * u - 2.0 * Ajφ * φ)

  return x
end

function get_identity_block(U::FESpace)
  return LinearAlgebra.I(num_free_dofs(U))
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
  solver=:block_cg,
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
Ωf = _fluid_mesh(model,params[:fluid][:domain])
V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
V_p = TestFESpace(Ωf,reffe_p;conformity=p_conformity(Ωf))
V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=BlockMultiFieldStyle())

# Trial spaces

z = zero(VectorValue{D,Float64})
u_bc = params[:bcs][:u][:values]
j_bc = params[:bcs][:j][:values]
U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
U_p = TrialFESpace(V_p)
U_φ = TrialFESpace(V_φ)
U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=BlockMultiFieldStyle())

# Weak form

res, jac = weak_form(params,k)
Tm = params[:matrix_type]
Tv = params[:vector_type]
assem = SparseMatrixAssembler(Tm,Tv,U,V)
op = FEOperator(res,jac,U,V,assem)

# Preconditioner

dΩj = Measure(get_triangulation(U_j),2*(k-1)+1)
Dj = assemble_matrix((j,v_j) -> ∫(j⋅v_j + (∇⋅j)⋅(∇⋅v_j))*dΩj ,U_j,V_j)

dΩp = Measure(get_triangulation(U_p),k-1)
Δp  = assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩp ,U_p,V_p)

Dj_solver = BackslashSolver()
Fu_solver = BackslashSolver()
Δp_solver = BackslashSolver()

P = MHDBlockPreconditioner(Dj_solver,Fu_solver,Δp_solver,Dj,Δp,params)
sysmat_solver = GMRESSolver(30,P,1e-6)

# Solve

xh = zero(U)
sysvec = residual(op,xh)
sysmat = jacobian(op,xh)

Pns = numerical_setup(symbolic_setup(P,sysmat),sysmat)
numerical_setup!(Pns,sysmat)

sysmat_ns = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)
numerical_setup!(sysmat_ns,sysmat)

x  = allocate_col_vector(sysmat)
dx = allocate_col_vector(sysmat)

solve!(x,sysmat_ns,sysvec)

abstol = 1e-6
for nliter in 1:10
  println("  > Newton Iteration: ", nliter)

  # Global linear solver
  println("    > Update linear solver: ")
  sysmat_solver = GMRESSolver(30,P,1e-6)
  sysmat_ns     = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)

  # Solve linearized problem
  println("    > Solver linearized problem:")
  rmul!(sysvec,-1)
  solve!(dx,sysmat_ns,sysvec)
  x .+= dx

  # Check convergence for the current residual
  xh = FEFunction(X,x)
  assemble_vector!(v->res(xh,v),sysvec,assem,Y)
  res_norm = norm(sysvec)
  (res_norm < abstol) && break;
  println("    > Residual: ", res_norm)

  # Assemble jacobian (fast in-place version)
  # and prepare solver
  assemble_matrix!((u,v)->jac(xh,u,v),sysmat,assem,X,Y)

end

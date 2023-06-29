using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra

using LinearAlgebra, BlockArrays

using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector

using GridapMHD: _hunt, add_default_params, _fluid_mesh, weak_form, _find_backend, p_conformity, _interior

function Gridap.Algebra.solve!(x::BlockVector,ns::GridapSolvers.LinearSolvers.GMRESNumericalSetup,b::BlockVector)
  solver, A, Pl, caches = ns.solver, ns.A, ns.Pl_ns, ns.caches
  m, tol = solver.m, solver.tol
  w, V, Z, H, g, c, s = caches
  println(" > Starting GMRES solver: ")

  # Initial residual
  mul!(w,A,x); w .= b .- w

  β    = norm(w)
  iter = 0
  while (β > tol)
    println("   > Iteration ", iter," - Residual: ", β)
    fill!(H,0.0)
    
    # Arnoldi process
    fill!(g,0.0); g[1] = β
    V[1] .= w ./ β
    j = 1
    while ( j < m+1 && β > tol )
      println("      > Inner iteration ", j," - Residual: ", β)
      # Arnoldi orthogonalization by Modified Gram-Schmidt
      solve!(Z[j],Pl,V[j])
      mul!(w,A,Z[j])
      for i in 1:j
        H[i,j] = dot(w,V[i])
        w .= w .- H[i,j] .* V[i]
      end
      H[j+1,j] = norm(w)
      V[j+1] = w ./ H[j+1,j]

      # Update QR
      for i in 1:j-1
        γ = c[i]*H[i,j] + s[i]*H[i+1,j]
        H[i+1,j] = -s[i]*H[i,j] + c[i]*H[i+1,j]
        H[i,j] = γ
      end

      # New Givens rotation, update QR and residual
      c[j], s[j], _ = LinearAlgebra.givensAlgorithm(H[j,j],H[j+1,j])
      H[j,j] = c[j]*H[j,j] + s[j]*H[j+1,j]; H[j+1,j] = 0.0
      g[j+1] = -s[j]*g[j]; g[j] = c[j]*g[j]

      β  = abs(g[j+1])
      j += 1
    end
    j = j-1

    # Solve least squares problem Hy = g by backward substitution
    for i in j:-1:1
      g[i] = (g[i] - dot(H[i,i+1:j],g[i+1:j])) / H[i,i]
    end

    # Update solution & residual
    for i in 1:j
      x .+= g[i] .* Z[i]
    end
    mul!(w,A,x); w .= b .- w
    println("        > Block residuals: ", map(norm,blocks(w)))

    iter += 1
  end
  println("   Exiting GMRES solver.")
  println("   > Num Iter: ", iter-1," - Final residual: ", β)

  return x
end

function Gridap.Algebra._check_convergence(nls,b,m0)
  m = Gridap.Algebra._inf_norm(b)
  println(">>>>>>>>>>>>>>>>>>>> Nonlinear Abs Error = $m")
  m < nls.tol * m0
end

"""
  This preconditioner is based on [(Li,2019)](https://doi.org/10.1137/19M1260372)
"""
struct MHDBlockPreconditioner <: Gridap.Algebra.LinearSolver
  Ij_solver
  Fu_solver
  Ip_solver
  Δp_solver
  Δφ_solver
  Ij
  Ip
  Δp
  Δφ
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
  Ij_ns
  Fu_ns
  Ip_ns
  Δp_ns
  Δφ_ns
  sysmat
  caches
end

function allocate_preconditioner_caches(solver::MHDBlockPreconditioner,A::BlockMatrix)
  du = allocate_col_vector(A[Block(1,1)])
  dp = allocate_col_vector(A[Block(2,2)])
  dj = allocate_col_vector(A[Block(3,3)])
  dφ = allocate_col_vector(A[Block(4,4)])
  return du, dp, dj, dφ
end

function Gridap.Algebra.numerical_setup(ss::MHDBlockPreconditionerSS, A::BlockMatrix)
  solver = ss.solver
  Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  Fk = Fu - (1.0/κ^2) * Kᵗ * K

  Ij_ns = numerical_setup(symbolic_setup(solver.Ij_solver,solver.Ij),solver.Ij)
  Fu_ns = numerical_setup(symbolic_setup(solver.Fu_solver,Fk),Fk)
  Ip_ns = numerical_setup(symbolic_setup(solver.Ip_solver,solver.Ip),solver.Ip)
  Δp_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Δp),solver.Δp)
  Δφ_ns = numerical_setup(symbolic_setup(solver.Δφ_solver,solver.Δφ),solver.Δφ)
  caches = allocate_preconditioner_caches(solver,A)
  return MHDBlockPreconditionerNS(ss.solver,Ij_ns,Fu_ns,Ip_ns,Δp_ns,Δφ_ns,A,caches)
end

function Gridap.Algebra.numerical_setup!(ns::MHDBlockPreconditionerNS, A::BlockMatrix)
  Fu = A[Block(1,1)]
  numerical_setup!(ns.Fu_ns,Fu)
  ns.sysmat = A
  return ns
end

# Follows Algorithm 4.1 in (Li,2019)
function Gridap.Algebra.solve!(x::BlockVector,ns::MHDBlockPreconditionerNS,b::BlockVector)
  caches = ns.caches
  bu, bp, bj, bφ = blocks(b)
  u, p, j, φ = blocks(x)
  du, dp, dj, dφ = caches

  # Solve for p
  #solve!(dp,ns.Ip_ns,bp)
  solve!(p,ns.Δp_ns,bp)
  #p .+= dp

  # Solve for φ
  solve!(φ,ns.Δφ_ns,bφ)

  # Solve for u
  solve!(u,ns.Fu_ns,bu)

  # Solve for j
  copy!(dj,bj)
  mul!(dj,sysmat[Block(3,1)],u,-1.0,1.0)
  solve!(j,ns.Ij_ns,dj)

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
params[:ζ] = nothing
res, jac = weak_form(params,k)
Tm = params[:matrix_type]
Tv = params[:vector_type]
assem = SparseMatrixAssembler(Tm,Tv,U,V)
op    = FEOperator(res,jac,U,V,assem)
al_op = Gridap.FESpaces.get_algebraic_operator(op)

# Preconditioner

Ω = Triangulation(model)
Γ = Boundary(model)
Λ = Skeleton(model)

dΩ = Measure(Ω,2*k)
dΓ = Measure(Γ,2*k)
dΛ = Measure(Λ,2*k)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

h_e_Λ = CellField(get_array(∫(1)dΛ),Λ)
h_e_Γ = CellField(get_array(∫(1)dΓ),Γ)

iRe = params[:fluid][:β]
γ = params[:fluid][:γ]
ζ = 100.0
aΛ(u,v) = ∫(-jump(u⋅n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ))*dΛ + ∫(ζ/h_e_Λ*jump(u⋅n_Λ)⋅jump(v⋅n_Λ))*dΛ
aΓ(u,v) = ∫(-(∇(u)⋅n_Γ)⋅v - u⋅(∇(v)⋅n_Γ))*dΓ + ∫(ζ/h_e_Γ*(u⋅n_Γ)⋅(v⋅n_Γ))*dΓ

ap(p,v_p) = (-γ)*(∫(∇(p)⋅∇(v_p))*dΩ + aΛ(p,v_p) + aΓ(p,v_p))
aφ(φ,v_φ) = (-γ)*(∫(∇(φ)⋅∇(v_φ))*dΩ + aΛ(φ,v_φ) + aΓ(φ,v_φ))

Ij = assemble_matrix((j,v_j) -> ∫(γ*j⋅v_j)*dΩ ,U_j,V_j)
Ip = assemble_matrix((p,v_p) -> ∫(iRe*p⋅v_p)*dΩ ,U_p,V_p)
Δp = assemble_matrix(ap,U_p,V_p)
Δφ = assemble_matrix(aφ,U_φ,V_φ)

Ij_solver = LUSolver()
Fu_solver = LUSolver()
Ip_solver = LUSolver()
Δp_solver = LUSolver()
Δφ_solver = LUSolver()

block_solvers = [Ij_solver,Fu_solver,Ip_solver,Δp_solver,Δφ_solver]
block_mats = [Ij,Ip,Δp,Δφ]
P = MHDBlockPreconditioner(block_solvers...,block_mats...,params)

sysmat_solver = GMRESSolver(500,P,1e-10)

# Gridap's Newton-Raphson solver
xh = zero(U)
sysvec = residual(op,xh)
sysmat = jacobian(op,xh)
sysmat_ns = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)

x  = allocate_col_vector(sysmat)
dx = allocate_col_vector(sysmat)

solve!(x,sysmat_ns,sysvec)

nlsolver = NewtonRaphsonSolver(sysmat_solver,1e-6,100)
nlsolver_cache = Gridap.Algebra.NewtonRaphsonCache(sysmat,sysvec,dx,sysmat_ns)
solve!(x,nlsolver,al_op,nlsolver_cache)

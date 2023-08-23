using FileIO
using BSON
using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra
using LinearAlgebra, BlockArrays
using GridapSolvers.LinearSolvers: allocate_col_vector, allocate_row_vector
using PartitionedArrays

using GridapMHD: weak_form, add_default_params


######################################################################################
# Preconditioner setup (same as BlockSolverTests.jl), cavity problem below

function Gridap.Algebra._check_convergence(nls,b,m0)
  m = Gridap.Algebra._inf_norm(b)
  println(">>>>>>>>>>>>>>>>>>>> Nonlinear Abs Error = $m")
  m < nls.tol * m0
end
function Gridap.Algebra._check_convergence(nls,b)
  m0 = Gridap.Algebra._inf_norm(b)
  println(">>>>>>>>>>>>>>>>>>>> Starting nonlinear solver")
  println(">>>>>>>>>>>>>>>>>>>> Nonlinear Abs Error = $m0")
  (false, m0)  
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
  Fk = Fu - (1.0/κ^2) * Kᵗ * K

  Dj_ns = numerical_setup(symbolic_setup(solver.Dj_solver,solver.Dj),solver.Dj)
  Fk_ns = numerical_setup(symbolic_setup(solver.Fk_solver,Fk),Fk)
  Δp_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Δp),solver.Δp)
  Ip_ns = numerical_setup(symbolic_setup(solver.Δp_solver,solver.Ip),solver.Ip)
  Iφ_ns = numerical_setup(symbolic_setup(solver.Iφ_solver,solver.Iφ),solver.Iφ)
  caches = allocate_preconditioner_caches(solver,A)
  return MHDBlockPreconditionerNS(ss.solver,Dj_ns,Fk_ns,Δp_ns,Ip_ns,Iφ_ns,A,caches)
end

function Gridap.Algebra.numerical_setup!(ns::MHDBlockPreconditionerNS, A::BlockMatrix)
  solver = ns.solver

  #! Pattern of matrix changes, so we need to recompute everything.
  # This will get fixed when we are using iterative solvers for Fk
  Fu = A[Block(1,1)]; K = A[Block(3,1)]; Kᵗ = A[Block(1,3)]; κ = solver.params[:fluid][:γ]
  Fk = Fu - (1.0/κ^2) * Kᵗ * K
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
  solve!(p,ns.Δp_ns,bp)
  solve!(dp,ns.Ip_ns,bp)
  p .= -α1 .* dp .- p

  #  Solve for φ
  #dφ .= -bφ
  solve!(φ,ns.Iφ_ns,bφ)

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

function BlockCavity(n; title="cavity", vtk=true)

  t_parts = get_part_ids(SequentialBackend(), 1)
  t = PTimer(t_parts, verbose=true)
  tic!(t, barrier=true)

  # Parameters
  ν = 1.0
  ρ = 1.0
  σ = 1.0
  B = VectorValue(0.0, 0.0, 10.0)
  f = VectorValue(0.0, 0.0, 0.0)
  u0 = 1.0
  B0 = norm(B)
  L = 1.0

  # Reduced quantities
  Re = u0 * L / ν
  Ha = B0 * L * sqrt(σ / (ρ * ν))
  N = Ha^2 / Re
  f̄ = (L / (ρ * u0^2)) * f
  B̄ = (1 / B0) * B
  α = 1.0
  β = 1.0 / Re
  γ = N

  # Domain and model
  #n=8
  domain = (0, L, 0, L, 0, L)
  cells = (n, n, n)
  model = simplexify(CartesianDiscreteModel(domain, cells))
  Ω = Interior(model)

  # Boundary conditions
  labels = get_face_labeling(model)
  Γw = append!(collect(1:4), [9, 10, 13, 14], collect(17:21), collect(23:26))
  Γl = append!(collect(5:8), [11, 12, 15, 16, 22])
  add_tag_from_tags!(labels, "wall", Γw)
  add_tag_from_tags!(labels, "lid", Γl)
  add_tag_from_tags!(labels, "insulating", "boundary")
  uw = VectorValue(0.0, 0.0, 0.0)
  ul = VectorValue(1.0, 0.0, 0.0)
  ji = VectorValue(0.0, 0.0, 0.0)

  _params = Dict(
      :ptimer => t,
      :debug => false,
      :solve => true,
      :res_assemble => false,
      :jac_assemble => false,
      :model => model,
      :fluid => Dict(
          :domain => model,
          :α => α,
          :β => β,
          :γ => γ,
          :f => f̄,
          :B => B̄,
      ),
      :bcs => Dict(
          :u => Dict(:tags => ["wall", "lid"], :values => [uw, ul]),
          :j => Dict(:tags => "insulating", :values => ji),
      )
  )

  params = add_default_params(_params)
  params[:ζ] = 1.0  # Augmented-Lagragian term 

  toc!(t, "pre_process")

  tic!(t; barrier=true)
  # ReferenceFEs
  k = 2
  T = Float64
  model = params[:model]
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian, VectorValue{D,T}, k)
  reffe_p = ReferenceFE(lagrangian, T, k - 1)
  reffe_j = ReferenceFE(raviart_thomas, T, k - 2)
  reffe_φ = ReferenceFE(lagrangian, T, k - 2)

  # Test spaces
  V_u = TestFESpace(model, reffe_u; dirichlet_tags=["wall", "lid"])
  V_p = TestFESpace(model, reffe_p; constraint=:zeromean)
  V_j = TestFESpace(model, reffe_j; dirichlet_tags="insulating")
  V_φ = TestFESpace(model, reffe_φ; conformity=:L2)
  #V = MultiFieldFESpace([V_u, V_p, V_j, V_φ])
  V  = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=BlockMultiFieldStyle())

  # Trial spaces
  U_u = TrialFESpace(V_u, [uw, ul])
  U_j = TrialFESpace(V_j, ji)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)
  #U = MultiFieldFESpace([U_u, U_p, U_j, U_φ])
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=BlockMultiFieldStyle())
  toc!(t, "fe_spaces")

  tic!(t; barrier=true)
  res, jac = weak_form(params, k)
  Tm = params[:matrix_type]
  Tv = params[:vector_type]
  assem = SparseMatrixAssembler(Tm, Tv, U, V)
  op = FEOperator(res, jac, U, V, assem)
  al_op = Gridap.FESpaces.get_algebraic_operator(op)

  dΩ = Measure(Ω,2*k)
  Dj = assemble_matrix((j,v_j) -> ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ ,U_j,V_j)
  Ij = assemble_matrix((j,v_j) -> ∫(j⋅v_j)*dΩ ,U_j,V_j)
  Δp = assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩ ,U_p,V_p)
  Ip = assemble_matrix((p,v_p) -> ∫(p*v_p)*dΩ,V_p,V_p)
  Iφ = assemble_matrix((φ,v_φ) -> ∫(-γ*φ*v_φ)*dΩ ,U_φ,V_φ)

  Dj_solver = LUSolver()
  Fk_solver = LUSolver()
  Δp_solver = LUSolver()
  Ip_solver = LUSolver()
  Iφ_solver = LUSolver()

  block_solvers = [Dj_solver,Fk_solver,Δp_solver,Ip_solver,Iφ_solver]
  block_mats = [Dj,Δp,Ip,Ij,Iφ]
  P = MHDBlockPreconditioner(block_solvers...,block_mats...,params)
  
  sysmat_solver = GMRESSolver(300,P,1e-8)

  # Gridap's Newton-Raphson solver
  xh = zero(U)
  sysmat = jacobian(op,xh)
  sysmat_ns = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)

  x  = allocate_col_vector(sysmat)
  dx = allocate_col_vector(sysmat)
  b  = allocate_col_vector(sysmat)

  A = allocate_jacobian(al_op,x)
  nlsolver = NewtonRaphsonSolver(sysmat_solver,1e-5,10)
  nlsolver_cache = Gridap.Algebra.NewtonRaphsonCache(A,b,dx,sysmat_ns)
  solve!(x,nlsolver,al_op,nlsolver_cache)

  toc!(t, "solve")

  if vtk
      tic!(t, barrier=true)
      ūh=FEFunction(U_u,x.blocks[1])
      p̄h=FEFunction(U_p,x.blocks[2])
      j̄h=FEFunction(U_j,x.blocks[3])
      φ̄h=FEFunction(U_φ,x.blocks[4])
      uh = u0 * ūh
      ph = (ρ * u0^2) * p̄h
      jh = (σ * u0 * B0) * j̄h
      φh = (u0 * B0 * L) * φ̄h
      writevtk(Ω, title, order=2, cellfields=["uh" => uh, "ph" => ph, "jh" => jh, "phi" => φh])
      toc!(t, "vtk")
  end

  info = Dict{Symbol,Any}()
  info[:ncells] = num_cells(model)
  info[:n_free_dofs_u] = length(x.blocks[1])
  info[:n_free_dofs_p] = length(x.blocks[2])
  info[:n_free_dofs_j] = length(x.blocks[3])
  info[:n_free_dofs_φ] = length(x.blocks[4])
  info[:n_dir_dofs_u] = length(get_dirichlet_dof_values(U_u))
  info[:n_dir_dofs_p] = length(get_dirichlet_dof_values(U_p))
  info[:n_dir_dofs_j] = length(get_dirichlet_dof_values(U_j))
  info[:n_dir_dofs_φ] = length(get_dirichlet_dof_values(U_φ))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  save("$title.bson", info)

end
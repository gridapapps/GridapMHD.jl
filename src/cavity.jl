
function cavity(;
  backend = nothing,
  np      = 1,
  title   = "Cavity",
  path    = ".",
  kwargs...)

  if isa(backend,Nothing)
    @assert np == 1
    info, t = _cavity(;title=title,path=path,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    if backend == :sequential
      info,t = with_debug() do distribute
        _cavity(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _cavity(;distribute=distribute,np=np,title=title,path=path,kwargs...)
      end
    end
  end
  info[:np] = np
  info[:backend] = backend
  info[:title] = title
  map_main(t.data) do data
    for (k,v) in data
      info[Symbol("time_$k")] = v.max
    end
    save(joinpath(path,"$title.bson"),info)
  end

  nothing
end

function _cavity(;
  distribute=nothing,
  np=1,
  nc=(4,4,4),
  ν=1.0,
  ρ=1.0,
  σ=1.0,
  B=VectorValue(0.0, 0.0, 10.0),
  f=VectorValue(0.0, 0.0, 0.0),
  L=1.0,
  u0=1.0,
  B0=norm(B),
  vtk=true,
  title="Cavity",
  path='.',
  solver=:julia,
  solver_params=default_solver_params(Val(solver)),
  verbose=true,
  )
  info = Dict{Symbol,Any}()

  @assert length(nc) == 3
  is_serial = isa(distribute,Nothing)

  if is_serial
    distribute = DebugArray
  end
  parts = distribute(LinearIndices((np,)))
  
  t = PTimer(parts,verbose=verbose)
  tic!(t, barrier=true)

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
  domain = (0, L, 0, L, 0, L)
  if is_serial
    model = simplexify(CartesianDiscreteModel(domain, nc))
  else
    # simplexify() is not implemented for DistributedDiscreteModel
    rank_partition = (np,1,1)
    model = CartesianDiscreteModel(parts,rank_partition,domain, nc)
  end
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
    :check_valid => false,
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
    ),
    :k => 2,
    :ζ => 0.0, # Augmented-Lagragian term
    :solver => solver,
    :solver_params => solver_params,
  )

  params = add_default_params(_params)
  toc!(t, "pre_process")

  tic!(t; barrier=true)
  # ReferenceFEs
  k = params[:k]
  T = Float64
  model = params[:model]
  D = num_cell_dims(model)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1)
  reffe_j = ReferenceFE(raviart_thomas,T,k-2)
  reffe_φ = ReferenceFE(lagrangian,T,k-2)

  mfs = (solver === :block_gmres) ? BlockMultiFieldStyle() : ConsecutiveMultiFieldStyle()

  # Test spaces
  V_u = TestFESpace(model, reffe_u; dirichlet_tags=["wall", "lid"])
  V_p = TestFESpace(model, reffe_p; constraint=:zeromean)
  V_j = TestFESpace(model, reffe_j; dirichlet_tags="insulating")
  V_φ = TestFESpace(model, reffe_φ; conformity=:L2)
  V   = MultiFieldFESpace([V_u, V_p, V_j, V_φ];style=mfs)

  # Trial spaces
  U_u = TrialFESpace(V_u, [uw, ul])
  U_j = TrialFESpace(V_j, ji)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)
  U   = MultiFieldFESpace([U_u, U_p, U_j, U_φ];style=mfs)
  toc!(t, "fe_spaces")

  tic!(t; barrier=true)
  res, jac = weak_form(params, k)
  Tm = params[:solver_params][:matrix_type]
  Tv = params[:solver_params][:vector_type]
  assem = SparseMatrixAssembler(Tm, Tv, U, V)
  op = FEOperator(res, jac, U, V, assem)

  if solver === :julia
    xh = zero(U)
    solver = NLSolver(show_trace=true, method=:newton)
    xh, solver_cache = solve!(xh, solver, op)
  elseif solver === :petsc
    petsc_options = params[:solver_params][:petsc_options]
    xh = GridapPETSc.with(args=split(petsc_options)) do
      xh = zero(U)
      solver = PETScNonlinearSolver()
      xh, solver_cache = solve!(xh, solver, op)
      snes_postpro(solver_cache,info)
      return xh
    end
  else
    xh = block_gmres_solver(parts,op,U,V,Ω,params)
  end
  toc!(t, "solve")

  if vtk
    tic!(t, barrier=true)
    ūh, p̄h, j̄h, φ̄h = xh
    uh = u0 * ūh
    ph = (ρ * u0^2) * p̄h
    jh = (σ * u0 * B0) * j̄h
    φh = (u0 * B0 * L) * φ̄h
    writevtk(Ω, joinpath(path,title), order=2, cellfields=["uh" => uh, "ph" => ph, "jh" => jh, "phi" => φh])
    toc!(t, "vtk")
  end

  info[:ncells]  = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(ūh))
  info[:ndofs_p] = length(get_free_dof_values(p̄h))
  info[:ndofs_j] = length(get_free_dof_values(j̄h))
  info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
  info[:ndofs]   = sum([info[:ndofs_u], info[:ndofs_p], info[:ndofs_j], info[:ndofs_φ]])
  info[:Re]      = Re
  info[:Ha]      = Ha

  return info, t
end

function block_gmres_solver(parts,op,U,V,Ω,params)
  petsc_solvers = [:mumps,:gmres_swartz,:amg,:cg_jacobi]
  uses_petsc = any(s -> s ∈ petsc_solvers, params[:solver_params][:block_solvers])
  if uses_petsc
    petsc_options = params[:solver_params][:petsc_options]
    xh = GridapPETSc.with(args=split(petsc_options)) do
      _block_gmres_solver(parts,op,U,V,Ω,params)
    end
  else
    xh = _block_gmres_solver(parts,op,U,V,Ω,params)
  end
  return xh
end

function _block_gmres_solver(parts,op,U,V,Ω,params)
  U_u, U_p, U_j, U_φ = U
  V_u, V_p, V_j, V_φ = V

  al_op = Gridap.FESpaces.get_algebraic_operator(op)
  xh = zero(U)
  sysmat = jacobian(op,xh)

  k  = params[:k]
  γ  = params[:fluid][:γ]
  dΩ = Measure(Ω,2*k)
  Dj = assemble_matrix((j,v_j) -> ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ ,U_j,V_j)
  Ij = assemble_matrix((j,v_j) -> ∫(j⋅v_j)*dΩ ,U_j,V_j)
  Δp = assemble_matrix((p,v_p) -> ∫(∇(p)⋅∇(v_p))*dΩ ,U_p,V_p)
  Ip = assemble_matrix((p,v_p) -> ∫(p*v_p)*dΩ,V_p,V_p)
  Iφ = assemble_matrix((φ,v_φ) -> ∫(-γ*φ*v_φ)*dΩ ,U_φ,V_φ)

  solver_params = params[:solver_params]
  block_solvers = map((s,m) -> get_block_solver(Val(s),m),solver_params[:block_solvers],[Dj,sysmat[Block(1,1)],Δp,Ip,Iφ])
  block_mats    = [Dj,Δp,Ip,Ij,Iφ]
  P = LI2019_Solver(block_solvers...,block_mats...,params)
  sysmat_solver = GMRESSolver(150,P,1e-8)

  test_block_solvers(parts,block_solvers,[Dj,sysmat[Block(1,1)],Δp,Ip,Iφ])

  # Gridap's Newton-Raphson solver
  sysmat_ns = numerical_setup(symbolic_setup(sysmat_solver,sysmat),sysmat)

  x  = allocate_col_vector(sysmat)
  dx = allocate_col_vector(sysmat)
  b  = allocate_col_vector(sysmat)

  A = allocate_jacobian(al_op,x)
  nlsolver = NewtonRaphsonSolver(sysmat_solver,1e-5,10)
  nlsolver_cache = Gridap.Algebra.NewtonRaphsonCache(A,b,dx,sysmat_ns)
  solve!(x,nlsolver,al_op,nlsolver_cache)

  ūh = FEFunction(U_u,x.blocks[1])
  p̄h = FEFunction(U_p,x.blocks[2])
  j̄h = FEFunction(U_j,x.blocks[3])
  φ̄h = FEFunction(U_φ,x.blocks[4])
  return [ūh, p̄h, j̄h, φ̄h]
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

function default_solver_params(::Val{:julia})
  return Dict(
    :matrix_type => SparseMatrixCSC{Float64,Int64},
    :vector_type => Vector{Float64},
  )
end

function default_solver_params(::Val{:petsc})
  Dict(
    :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type   => Vector{PetscScalar},
    :petsc_options => "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  )
end

function default_solver_params(::Val{:block_gmres})
  Dict(
    :matrix_type   => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type   => Vector{PetscScalar},
    :block_solvers => [:mumps,:mumps,:mumps,:mumps,:mumps],
    :petsc_options => "-ksp_error_if_not_converged true -ksp_converged_reason"
  )
end

get_block_solver(::Val{:julia},m) = LUSolver()
get_block_solver(::Val{:mumps},m) = PETScLinearSolver(cavity_mumps_setup)
get_block_solver(::Val{:amg},m) = PETScLinearSolver(cavity_amg_setup)
get_block_solver(::Val{:gmres_swartz},m) = PETScLinearSolver(cavity_gmres_setup)
get_block_solver(::Val{:cg_jacobi},m) = PETScLinearSolver(cavity_cg_setup)

#function get_block_solver(::Val{:cg_jacobi},m)
#  P = RichardsonSmoother(JacobiLinearSolver(),10,2.0/3.0)
#  Pns = numerical_setup(symbolic_setup(P,m),m)
#  return IS_ConjugateGradientSolver(;Pl=Pns,maxiter=10)
#end

function cavity_mumps_setup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
  # percentage increase in the estimated working space
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  14, 1000)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-6)
end

# Two iterations of standalone AMG solver
function cavity_amg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(2)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  #@check_error_code GridapPETSc.PETSC.PCGAMGSetType(pc[],GridapPETSc.PETSC.PCGAMGAGG) # or PCGAMGCLASSICAL
  #@check_error_code GridapPETSc.PETSC.PCGAMGSetNSmooths(pc[],0)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
end

# GMRES + Additive Schwartz preconditioner
function cavity_gmres_setup(ksp)
  rtol = PetscScalar(1.e-3)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  # GMRES solver
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)

  # Additive Schwartz preconditioner
  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCASM)
  #@check_error_code GridapPETSc.PETSC.PCASMSetLocalType(pc[],GridapPETSc.PETSC.PC_COMPOSITE_ADDITIVE)
end

# CG + Jacobi preconditioner, 10 iterations
function cavity_cg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(10)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPCG)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCJACOBI)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
end

############################################################################################
# FIXES

# Overload that will be removed when bug is fixed in PartitionedArrays
function Base.similar(a::PSparseMatrix,::Type{T},inds::Tuple{<:PRange,<:PRange}) where T
  rows,cols = inds
  matrix_partition = map(partition(a),partition(rows),partition(cols)) do values, row_indices, col_indices
    PartitionedArrays.allocate_local_values(values,T,row_indices,col_indices)
  end
  PSparseMatrix(matrix_partition,partition(rows),partition(cols))
end

function Base.similar(::Type{<:PSparseMatrix{V}},inds::Tuple{<:PRange,<:PRange}) where V
  rows,cols = inds
  matrix_partition = map(partition(rows),partition(cols)) do row_indices, col_indices
    PartitionedArrays.allocate_local_values(V,row_indices,col_indices)
  end
  PSparseMatrix(matrix_partition,partition(rows),partition(cols))
end

function Base.copy(a::PSparseMatrix)
  matrix_partition = similar(a.matrix_partition)
  copy!(matrix_partition, a.matrix_partition)
  PSparseMatrix(matrix_partition,a.row_partition,a.col_partition)
end

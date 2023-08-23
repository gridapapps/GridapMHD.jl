
function cavity(;
  backend = nothing,
  np      = nothing,
  title   = "Cavity",
  path    = ".",
  kwargs...)

  if isa(backend,Nothing)
    @assert isa(np,Nothing)
    info, t = _cavity(;title=title,path=path,mesh=mesh,kwargs...)
  else
    @assert backend ∈ [:sequential,:mpi]
    @assert !isa(np,Nothing)
    np = isa(np,Int) ? (np,) : np
    if backend == :sequential
      info,t = with_debug() do distribute
        _cavity(;distribute=distribute,rank_partition=np,title=_title,path=path,mesh=mesh,kwargs...)
      end
    else
      info,t = with_mpi() do distribute
        _cavity(;distribute=distribute,rank_partition=np,title=_title,path=path,mesh=mesh,kwargs...)
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
  rank_partition=nothing,
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
  solver=:julia,
  )

  @assert length(nc) == 3
  is_serial = isa(distribute,Nothing)

  if is_serial
    @assert isa(rank_partition,Nothing)
    rank_partition = (1,)
    distribute = DebugArray
  end
  @assert length(rank_partition) == length(nc)
  parts = distribute(LinearIndices((prod(rank_partition),)))
  
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
    model = simplexify(CartesianDiscreteModel(parts,rank_partition,domain, nc))
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
      :ζ => 0.0 # Augmented-Lagragian term 
  )

  params = add_default_params(_params)
  toc!(t, "pre_process")

  tic!(t; barrier=true)
  # ReferenceFEs
  k = 2
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
  Tm = params[:matrix_type]
  Tv = params[:vector_type]
  assem = SparseMatrixAssembler(Tm, Tv, U, V)
  op = FEOperator(res, jac, U, V, assem)

  if solver === :block_gmres
    xh = block_gmres_solver(op,U,V)
  else
    xh = zero(U)
    solver = NLSolver(show_trace=true, method=:newton)
    xh, cache = solve!(xh, solver, op)
  end
  toc!(t, "solve")

  if vtk
    tic!(t, barrier=true)
    ūh, p̄h, j̄h, φ̄h = xh
    uh = u0 * ūh
    ph = (ρ * u0^2) * p̄h
    jh = (σ * u0 * B0) * j̄h
    φh = (u0 * B0 * L) * φ̄h
    writevtk(Ω, title, order=2, cellfields=["uh" => uh, "ph" => ph, "jh" => jh, "phi" => φh])
    toc!(t, "vtk")
  end

  info = Dict{Symbol,Any}()
  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(ūh))
  info[:ndofs_p] = length(get_free_dof_values(p̄h))
  info[:ndofs_j] = length(get_free_dof_values(j̄h))
  info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  save("$title.bson", info)
end
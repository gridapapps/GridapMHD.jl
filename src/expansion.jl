function expansion(;
  backend=nothing,
  np=nothing,
  title = "Expansion",
  mesh = "710",
  path=".",
  kwargs...)

    if backend === nothing
      @assert np === nothing
      info, t = _expansion(;title=title,path=path,mesh=mesh,kwargs...)
    else
      @assert backend !== nothing
      @assert np !== nothing
      info, t = with_backend(_find_backend(backend),np) do _parts
        _expansion(;parts=_parts,title=title,path=path,mesh=mesh,kwargs...)
      end
    end
    info[:np] = np
    info[:backend] = backend
    info[:title] = title
    info[:mesh] = mesh
    map_main(t.data) do data
      for (k,v) in data
        info[Symbol("time_$k")] = v.max
      end
      save(joinpath(path,"$title.bson"),info)
    end

  nothing
end

function _expansion(;
  parts=nothing,
  title = "Expansion",
  mesh = "720",
  vtk=true,
  path=".",
  debug = false,
  verbose=true,
  solver=:julia,
  Z = 4.0, #Expansion Ratio, it has to be consistent with the mesh
  b = 0.2, #Outlet channel aspect ratio
  N = 1.0,
  Ha = 1.0,
  cw = 0.028,
  τ = 100,
  inlet = :parabolic,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0"
)
  
  info = Dict{Symbol,Any}()

  if parts === nothing
    t_parts = get_part_ids(SequentialBackend(),1)
  else
    t_parts = parts
  end
  t = PTimer(t_parts,verbose=verbose)
  tic!(t,barrier=true)

  # The domain is of size L_out x 2 x 2 and L_in x 2/Z x 2
  # after and before the expansion respectively.
  msh_file = joinpath(@__FILE__,"..","..","meshes","expansion","Expansion_"*mesh*".msh") |> normpath
  model = GmshDiscreteModel(parts,msh_file)
  if debug && vtk
    writevtk(model,"expansion_model")
    toc!(t,"model")
  end
  Ω = Interior(model,tags="PbLi")
  toc!(t,"triangulation")

  # Parameters and bounday conditions

  # Option 1 (CFD)
  # α = 1.0
  # β = 1.0/Re
  # γ = N
  # Option 2 (MHD) is chosen in the experimental article
  α = (1.0/N)
  β = (1.0/Ha^2)
  γ = 1.0
  
# This gives a parabolic profile with mean(u_inlet)=Z , which ensures mean(u_outlet) = 1
#  u_inlet_parabolic((x,y,z)) = VectorValue(36.0*Z*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0)
   
# Shercliff velocity profile at the entrance with mean(u_inlet) = Z, which ensures mean(u_outlet) = 1
#  kp_inlet = GridapMHD.kp_shercliff_cartesian(1/Z,Ha)
#  u_inlet_shercliff((x,y,z)) = VectorValue(GridapMHD.analytical_GeneralHunt_u(1.0/Z, 0.0, -Z*kp_inlet, Ha,200,(z,y,x))[3],0.0,0.0)
  
#  u_inlet_cte = VectorValue(Z,0.0,0.0)
  
  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :solve=>true,
    :res_assemble=>false,
    :jac_assemble=>false,
    :model => model,
    :fluid=>Dict(
      :domain=>model,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :f=>VectorValue(0.0,0.0,0.0),
      :B=>VectorValue(0.0,1.0,0.0),
    )
   )

  if cw == 0.0
   params[:bcs] = Dict( 
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet(inlet,Ha,Z,b), VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
		:tags => ["wall", "inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
      )
    )

  else 
   params[:bcs] = Dict(
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet(inlet,Ha,Z,b), VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
        :tags => ["inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
      ),
      :thin_wall => [Dict(
        :τ=>τ,
        :cw=>cw,
        :domain => Boundary(model, tags="wall")
      )]
    )
  end

  toc!(t,"pre_process")

  # Solve it
  if solver == :julia
    params[:solver] = NLSolver(show_trace=true,method=:newton)
    xh,fullparams = main(params)
  elseif solver == :petsc
    xh,fullparams = GridapPETSc.with(args=split(petsc_options)) do
    params[:matrix_type] = SparseMatrixCSR{0,PetscScalar,PetscInt}
    params[:vector_type] = Vector{PetscScalar}
    params[:solver] = PETScNonlinearSolver()
    params[:solver_postpro] = cache -> snes_postpro(cache,info)
    xh = main(params)
    end
  else
    error()
  end
  t = fullparams[:ptimer]
  tic!(t,barrier=true)

  uh,ph,jh,φh = xh
  div_jh = ∇·jh
  div_uh = ∇·uh
  Grad_p = ∇·ph

  if vtk
    writevtk(Ω,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,"div_uh"=>div_uh,"div_jh"=>div_jh,"kp"=>Grad_p])
    toc!(t,"vtk")
  end
  if verbose
    display(t)
  end

  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(uh))
  info[:ndofs_p] = length(get_free_dof_values(ph))
  info[:ndofs_j] = length(get_free_dof_values(jh))
  info[:ndofs_φ] = length(get_free_dof_values(φh))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Ha] = Ha
  info[:N]  = N
  info[:Re] = Ha^2/N
  info[:cw] = cw
  info[:τ] = τ
  info, t

end

function u_inlet(inlet,Ha,Z,b)
 
  u_inlet_parabolic((x,y,z)) = VectorValue(36.0*Z*(y-1/Z)*(y+1/Z)*(z-b*Z)*(z+b*Z),0,0)  
    
  kp_inlet = GridapMHD.kp_shercliff_cartesian(b*Z,Ha/Z)
  u_inlet_shercliff((x,y,z)) = VectorValue(GridapMHD.analytical_GeneralHunt_u(b*Z, 0.0, -Z*kp_inlet, Ha/Z,200,(z*Z,y*Z,x))[3],0.0,0.0) 

  u_inlet_cte = VectorValue(Z,0.0,0.0)  

  if inlet == :parabolic
	U = u_inlet_parabolic
  elseif inlet == :shercliff
        U = u_inlet_shercliff
  else
        U = u_inlet_cte
  end
 U
end

	

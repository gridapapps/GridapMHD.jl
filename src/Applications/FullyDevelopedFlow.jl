function FullyDeveloped(;
  backend=nothing,
  np=nothing,
  parts=nothing,
  title = "FullyDeveloped",
  nruns=1,
  path=".",
  kwargs...)

  for ir in 1:nruns
    _title = title*"_r$ir"
    @assert parts === nothing
    if backend === nothing
      @assert np === nothing
      info, t = _FullyDeveloped(;title=_title,path=path,kwargs...)
    else
      @assert backend !== nothing
      info, t = with_backend(_find_backend(backend),(np...,1)) do _parts
        _FullyDeveloped(;parts=_parts,title=_title,path=path,kwargs...)
      end
    end
    info[:np] = np
    info[:backend] = backend
    info[:title] = title
    map_main(t.data) do data
      for (k,v) in data
        info[Symbol("time_$k")] = v.max
      end
      save(joinpath(path,"$_title.bson"),info)
    end
  end

  nothing
end

function _FullyDeveloped(;
  parts=nothing,
  nc=(3,3),
  Ha = 10.0,
  dir_B=(0.0,1.0,0.0),
  b=1.0,
  vtk=true,
  title="test",
  path=".",
  debug=false,
  res_assemble=false,
  jac_assemble=false,
  solve=true,
  solver=:julia,
  verbose=true,
  mesh = false,
  strech_fine = false,
  cw_Ha = 0.0,
  cw_s = 0.0,
  τ_Ha = 100.0,
  τ_s = 100.0,
  nsums = 10,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  )

  info = Dict{Symbol,Any}()

  if parts === nothing
    t_parts = get_part_ids(SequentialBackend(),1)
  else
    t_parts = parts
  end
  t = PTimer(t_parts,verbose=verbose)
  tic!(t,barrier=true)

  domain = (-b,b,-1.0,1.0,0.0,0.1)

  # Reduced quantities
  Re = 1.0
  N = Ha^2/Re
  dirB = (1/norm(VectorValue(dir_B)))*VectorValue(dir_B)
  α = 1.0/N
  β = 1.0/Ha^2
  γ = 1.0

  # Prepare problem in terms of reduced quantities

   strech_Ha = sqrt(Ha/(Ha-1))
   strech_side = sqrt(sqrt(Ha)/(sqrt(Ha)-1))

  function map1(coord)
     ncoord = GridapMHD.strechMHD(coord,domain=(0,-b,0,-1.0),factor=(strech_side,strech_Ha),dirs=(1,2))
     ncoord = GridapMHD.strechMHD(ncoord,domain=(0,b,0,1.0),factor=(strech_side,strech_Ha),dirs=(1,2))
     ncoord  
   end

  function map_fine(coord)
     ncoord = GridapMHD.strechMHD(coord,domain=(0,-b,0,-1.0),factor=(strech_Ha,strech_Ha),dirs=(1,2))
     ncoord = GridapMHD.strechMHD(ncoord,domain=(0,b,0,1.0),factor=(strech_Ha,strech_Ha),dirs=(1,2))
     ncoord
  end

 partition=(nc[1],nc[2],3)
 if strech_fine
 model = CartesianDiscreteModel(
     parts,domain,partition;isperiodic=(false,false,true),map=map_fine)
 else
 model = CartesianDiscreteModel(
     parts,domain,partition;isperiodic=(false,false,true),map=map1)
 end
  Ω = Interior(model)
 
  labels = get_face_labeling(model)
  tags_j_Ha = append!(collect(1:20),[23,24])
  tags_j_side = append!(collect(1:20),[25,26])
  tags_outlet = append!(collect(1:20),[22])

  add_tag_from_tags!(labels,"Ha_walls",tags_j_Ha)
  add_tag_from_tags!(labels,"side_walls",tags_j_side)
  add_tag_from_tags!(labels,"outlet",tags_outlet)
  
  if mesh
    writevtk(model,"Mesh")
  end

  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :model=>model,
    :res_assemble=>res_assemble,
    :jac_assemble=>jac_assemble,
    :solve=>solve,
    :fluid=>Dict(
      :domain=>model,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :f=>VectorValue(0.0,0.0,1.0),
      :B=>dirB,
    ),
    :bcs => Dict(
      :u=>Dict(:tags=>["Ha_walls","side_walls"]),
      :j=>Dict{Symbol,Vector{String}}(),
      :thin_wall=>Vector{Dict{Symbol,Any}}() 
    )
  )

  insulated_tags = Vector{String}()
  thinWall_options = Vector{Dict{Symbol,Any}}() 

  if cw_Ha == 0.0 
    push!(insulated_tags,"Ha_walls")
  else
    push!(thinWall_options, Dict(
        :cw => cw_Ha,
        :τ => τ_Ha, 
        :domain =>Boundary(model, tags="Ha_walls")
        )
      ) 
  end  

  if cw_s == 0.0
    push!(insulated_tags,"side_walls")
  else
    push!(thinWall_options, Dict(
        :cw => cw_s,
        :τ => τ_s,
        :domain =>Boundary(model, tags="side_walls")
        )
      )
  end
  
#  if insulated_tags != [] 
    params[:bcs][:j] = Dict(:tags=>insulated_tags)  
#  end

  if thinWall_options != []
    params[:bcs][:thin_wall] = thinWall_options
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
      main(params)
    end
  else
    error()
  end
  t = fullparams[:ptimer]

  # Rescale quantities

  tic!(t,barrier=true)
  uh,ph,jh,φh = xh
  
  div_jh = ∇·jh 
  div_uh = ∇·uh  
  
  
  #Post process
  
  dΩ = Measure(Ω,6)
  uh_0 = sum(∫(uh)*dΩ)[3]/sum(∫(1.0)*dΩ)
 
  uh_n = uh/uh_0
  ph_n = ph/uh_0
  jh_n = jh/uh_0
  φh_n  =φh/uh_0  

  div_jh_n = div_jh/uh_0
  div_uh_n = div_uh/uh_0

  if cw_s == 0.0
    u_a(x) = analytical_GeneralHunt_u(1.0,cw_Ha,-1.0,Ha,nsums,x)
    e_u = u_a - uh_n
  else
    e_u =  uh_n
  end

  kp = 1/uh_0

  if cw_s == 0.0 && cw_Ha == 0.0
    kp_a = kp_shercliff_cartesian(b,Ha)
  else
    kp_a = kp_tillac(b,Ha,cw_s,cw_Ha)
  end

  dev_kp = 100*abs(kp_a-kp)/max(kp_a,kp)

  if vtk
    writevtk(Ω,joinpath(path,title),
      order=2,
      cellfields=[
        "uh"=>uh_n,"e_u"=>e_u,"ph"=>ph_n,"jh"=>jh_n,"phi"=>φh_n,"div_jh"=>div_jh_n,"div_uh"=>div_uh_n])
    toc!(t,"vtk")
  end
  if verbose
    display(t)
  end

  info[:nc] = nc
  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(uh))
  info[:ndofs_p] = length(get_free_dof_values(ph))
  info[:ndofs_j] = length(get_free_dof_values(jh))
  info[:ndofs_φ] = length(get_free_dof_values(φh))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  info[:N] = N
  info[:cw_s] = cw_s
  info[:τ_s] = τ_s
  info[:cw_Ha] = cw_Ha
  info[:τ_Ha] = τ_Ha
  info[:b] = b
  info[:uh_0] = uh_0
  info[:kp] = kp
  info[:kp_a] =kp_a
  info[:dev_kp] = dev_kp 
  info, t
end

function analytical_GeneralHunt_u(
  #General Hunt analytical formula (d_b = 0 for Shercliff)
  l::Float64,       # channel aspect ratio
  d_b::Float64,     # Hartmann walls conductivity ratio
  grad_pz::Float64, # Dimensionless (MHD version) presure gradient
  Ha::Float64,      # Hartmann number
  n::Int,           # number of sumands included in Fourier series
  x)                # evaluation point normaliced by the Hartmann characteristic lenght

  V = 0.0; V0=0.0;
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = (Ha^2 + 4*α_k^2)^(0.5)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)
    
    eplus_1k = 1 + exp(-2*r1_k)
    eminus_1k = 1 - exp(-2*r1_k)
    eplus_2k = 1 + exp(-2*r2_k)
    eminus_2k = 1 - exp(-2*r2_k)
    eplus_k = 1 + exp(-2*(r1_k+r2_k))
    e_x_1k = 0.5*(exp(-r1_k*(1-x[2]))+exp(-r1_k*(1+x[2])))
    e_x_2k = 0.5*(exp(-r2_k*(1-x[2]))+exp(-r2_k*(1+x[2])))
    
    V2 = ((d_b*r2_k + eminus_2k/eplus_2k)*e_x_1k)/(0.5*N*d_b*eplus_1k + eplus_k/eplus_2k)
    V3 = ((d_b*r1_k + eminus_1k/eplus_1k)*e_x_2k)/(0.5*N*d_b*eplus_2k + eplus_k/eplus_1k)

    V += 2*(-1)^k*cos(α_k * x[1])/(l*α_k^3) * (1-V2-V3)
  end
  u_z = V*Ha^2*(-grad_pz) 

  VectorValue(0.0*u_z,0.0*u_z,u_z)
end

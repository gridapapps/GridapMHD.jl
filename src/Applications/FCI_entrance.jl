# Change of the electrical insulation at the entrance of an Flow Channel Insert (FCI)
# Validation problem with respect to the experiment published by L. Bühler et al. "Experimental investigation of liquid metal MHD flow entering a flow channel insert"
# Fusion Engineering and Design 154(2020):111484

function FCI_Entrance(;
  backend=nothing,
  np=nothing,
  title = "FCI_Entrance",
  mesh = "Mesh",
  path = datadir(),
  kwargs...
  )

  if backend === nothing
    @assert np === nothing
    info, t = _FCI_Entrance(;title=title,path=path,mesh=mesh,kwargs...)
  else
    @assert backend !== nothing
    @assert np !== nothing
    info, t = prun(_find_backend(backend),np) do _parts
      _FCI_Entrance(;parts=_parts,title=title,path=path,mesh=mesh,kwargs...)
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



function _FCI_Entrance(;	
  Ha = 1.0,               # Hartmann number
  Re = 1.0,               # Reynolds number
  c_w = 0.0727,           # Wall conductance ratio 
  c_FCI = 0.00476,        # FCI conductance ratio
  B_dir = (1.0,1.0,0.0),  # Direction of the B-field (norm =1)  
  τ_w = 100,              # Penalty for the wall boundary condition
  τ_FCI =100,             # Penalty for the FCI boundary condition
  R = 1,                  # Radious of the pipe
  ρ = 1,                  # Density (Kg/m3)
  ν = 1,                  # Kinematic viscosity (Pa·s)
  σ = 1,                  # Electrical conductivity (S/m)
  dimensionless = true,   # Results with dimensionless variables
  distribute = nothing,   # For the parallel call
  rank_partition = nothing,   # Processor layout
  debug = false,          # debugging mode
  vtk = true,             # vtk file as an output
  mesh   = "Mesh",          # Mesh file .msh
  title  = "FCI_entrance", # Title of the job
  path   = datadir(),    # path for writting the results
  solver = :petsc,        # solver used for the analyses (julia or petsc)
  )
#------------- Derived variables------------------ 
  
  N = Ha^2/Re # Interaction parameter or Stuart number
  
  B_dir = (1/norm(VectorValue(B_dir)))*VectorValue(B_dir) #Direction of the magnetic field normalized

  U0 = Re*ν/R  # Velocity scale	
  B0 = Ha*sqrt((ρ*ν)/σ)/R 	   # Module of the external B-field	

#------------Mesh and geometry------------------------

  info = Dict{Symbol,Any}()
  
  if isa(distribute,Nothing)
    @assert isa(rank_partition,Nothing)
    parts = DebugArray(collect(LinearIndices((1,))))
  else
    parts = distribute(LinearIndices((prod(rank_partition),)))
  end

  t = PTimer(parts,verbose=true)
  tic!(t,barrier=true)

  msh_file = joinpath(@__FILE__,"..","..","meshes","FCI_"*mesh*".msh") |> normpath 	#find the FCI_*mesh".msh file in the "meshes" folder
  model = GmshDiscreteModel(parts,msh_file)
  
  if debug && vtk
    writevtk(model,"MeshDebug")
    toc!(t,"model")
  end 
  
  Ω = Triangulation(model,tags="Fluid")
   
  toc!(t,"triangulation")

#-----------Boundary Conditions-----------------------

  z = VectorValue(0.0, 0.0, 0.0)

  # Parabolic profile at the inlet
  function g_u(coord)
     r = sqrt(coord[1]^2 + coord[2]^2)
     return VectorValue(0.0, 0.0,2*U0*(R^2-r^2)/R^2)
    
  end
  
  #Thin wall boundaries

  Γ_w   = Boundary(model, tags="wall")
  Γ_FCI = Boundary(model, tags="FCI")

#-------Dictionary of parameters for the solver-----------

  if isa(solver,Symbol)
    solver = default_solver_params(Val(solver))
  end

  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :fluid=>Dict(
      :domain=>model,
      :α=>1/N,
      :β=>1/Ha^2,
      :γ=>1,
      :u=>Dict(
        :tags=>["wall","FCI","inlet"],
        :values=>[z,z,g_u]),
      :f=>z,
      :B=>B_dir,
    ),
    :solver=>solver,
  )

  if c_w == 0 && c_FCI == 0
    params[:fluid][:j] = Dict(
      :tags=>["inlet","outlet","FCI"],
      :values=>[z,z,z])
  elseif c_w != 0 && c_FCI == 0
    params[:fluid][:j] = Dict(
      :tags=>["inlet","outlet","FCI"],
      :values=>[z,z,z])
    params[:fluid][:thin_wall] = Dict(
      :domain=>Γ_w,
      :cw=>c_w,
      :τ=>τ_w,
      :jw=>0)
  elseif c_w == 0 && c_FCI != 0
    params[:fluid][:j]= Dict(
      :tags=>["inlet","outlet","wall"],
      :values=>[z,z,z])
    params[:fluid][:thin_wall] = Dict(
      :domain=>Γ_FCI,
      :cw=>c_FCI,
      :τ=>τ_FCI,
      :jw=>0)
  else
    params[:fluid][:j]= Dict(
      :tags=>["inlet","outlet"],
      :values=>[z,z])
    params[:fluid][:thin_wall] = [
      Dict(
        :domain=>Γ_FCI,
        :cw=>c_FCI,
        :τ=>τ_FCI,
        :jw=>0)
      ,
      Dict(
        :domain=>Γ_w,
        :cw=>c_w,
        :τ=>τ_w,
        :jw=>0)
    ]
  end

  toc!(t,"pre_process")

#-----------------Solver call-------------------------

  if !uses_petsc(Val(params[:solver][:solver]),params[:solver])
    xh,fullparams,info = main(params;output=info)
  else
    petsc_options = params[:solver][:petsc_options]
    xh,fullparams,info = GridapPETSc.with(args=split(petsc_options)) do
      main(params;output=info)
    end
  end

#---------------Post process-------------------------

  t = params[:ptimer]
  tic!(t,barrier=true)
  
  if dimensionless
   uh,ph,jh,φh = xh
  else
   ūh,p̄h,j̄h,φ̄h = xh
   uh = U0*ūh
   ph = (σ*U0*L*B0^2)*p̄h
   jh = (σ*U0*B0)*j̄h
   φh = (U0*B0*L)*φ̄h
  end

  divj = ∇⋅jh

  if vtk
    writevtk(Ω,joinpath(path,title),order=2,cellfields=["u"=>uh,"p"=>ph,"j"=>jh,"phi"=>φh,"divj"=>divj])
    toc!(t,"vtk")
  end

  info[:ncells] = num_cells(model)
  info[:ndofs_u] = length(get_free_dof_values(uh))
  info[:ndofs_p] = length(get_free_dof_values(ph))
  info[:ndofs_j] = length(get_free_dof_values(jh))
  info[:ndofs_φ] = length(get_free_dof_values(φh))
  info[:ndofs] =   length(get_free_dof_values(xh))
  info[:Ha] = Ha
  info[:N]  = N
  info[:Re] = Ha^2/N
  info, t
end



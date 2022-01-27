function hunt(;
  nx=3,
  ny=3,
  ν=1.0,
  ρ=1.0,
  σ=1.0,
  B=VectorValue(0.0,1.0,0.0),
  f=VectorValue(0.0,0.0,1.0),
  L=1.0,
  u0=1.0,
  B0=norm(B),
  nsums = 10,
  vtk=true,
  title="test",
  debug=false,
  solver="julia",
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
  )

  t = PTimer(get_part_ids(sequential,1),verbose=true)
  tic!(t,barrier=true)

  domain_phys = (-L,L,-L,L,0.0*L,0.1*L)

  # Reduced quantities
  Re = u0*L/ν
  Ha = B0*L*sqrt(σ/(ρ*ν))
  N = Ha^2/Re
  f̄ = (L/(ρ*u0^2))*f
  B̄ = (1/B0)*B
  α = 1.0
  β = 1.0/Re
  γ = N
  domain = domain_phys ./ L

  # Prepare problem in terms of reduced quantities
  map(x) = (2/sqrt(2))*VectorValue(
    sign(x[1])*(abs(x[1])*0.5)^0.5,
    sign(x[2])*(abs(x[2])*0.5)^0.5,
    x[3]*sqrt(2)/2)
  partition=(nx,ny,3)
  model = CartesianDiscreteModel(
    domain,partition;isperiodic=(false,false,true),map=map)
  Ω = Interior(model)
  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j = append!(collect(1:20),[25,26])
  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"insulating",tags_j)

  params = Dict(
    :ptimer=>t,
    :debug=>debug,
    :fluid=>Dict(
      :domain=>Ω,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :u=>Dict(
        :tags=>"noslip",
        :values=>VectorValue(0,0,0)),
      :j=>Dict(
        :tags=>"insulating",
        :values=>VectorValue(0,0,0)),
      :f=>f̄,
      :B=>B̄,
    ),
  )

  toc!(t,"pre_process")

  # Solve it
  if solver == "julia"
    params[:solver] = NLSolver(show_trace=true,method=:newton)
    xh = main(params)
  elseif solver == "petsc"
    xh = GridapPETSc.with(args=split(petsc_options)) do
    params[:matrix_type] = SparseMatrixCSR{0,PetscScalar,PetscInt}
    params[:vector_type] = Vector{PetscScalar}
    params[:solver] = PETScNonlinearSolver()
    xh = main(params)
    end
  else
    error()
  end
  t = params[:ptimer]

  # Rescale quantities

  tic!(t,barrier=true)
  ūh,p̄h,j̄h,φ̄h = xh
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (u0*B0*L)*φ̄h
  grid_phys = UnstructuredGrid(get_grid(model))
  node_coords = get_node_coordinates(grid_phys)
  node_coords .= L .* node_coords
  Ω_phys = Gridap.Geometry.BodyFittedTriangulation(Ω.model,grid_phys,Ω.tface_to_mface)

  # Post process

  μ = ρ*ν
  grad_pz = -f[3]/ρ
  u(x) = analytical_hunt_u(L,L,μ,grad_pz,Ha,nsums,x)
  j(x) = analytical_hunt_j(L,L,σ,μ,grad_pz,Ha,nsums,x)

  if vtk
    writevtk(Ω_phys,"$(title)_Ω_fluid",
      cellfields=[
      "uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,"u"=>u,"j"=>j,])
  end

  k = 2
  dΩ_phys = Measure(Ω_phys,2*k)
  eu = u - uh
  ej = j - jh
  eu_h1 = sqrt(sum(∫( ∇(eu)⊙∇(eu) + eu⋅eu  )dΩ_phys))
  eu_l2 = sqrt(sum(∫( eu⋅eu )dΩ_phys))
  ej_l2 = sqrt(sum(∫( ej⋅ej )dΩ_phys))
  toc!(t,"post_process")
  display(t)

  info = Dict{Symbol,Any}()
  info[:ncells_fluid] = num_cells(Ω)
  info[:ndofs_u] = length(get_free_dof_values(ūh))
  info[:ndofs_p] = length(get_free_dof_values(p̄h))
  info[:ndofs_j] = length(get_free_dof_values(j̄h))
  info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
  info[:ndofs] = length(get_free_dof_values(xh))
  info[:Re] = Re
  info[:Ha] = Ha
  info[:eu_h1] = eu_h1
  info[:eu_l2] = eu_l2
  info[:ej_l2] = ej_l2

  info
end

function analytical_hunt_u(
  a::Float64,       # semi-length of side walls
  b::Float64,       # semi-length of Hartmann walls
  μ::Float64,       # fluid viscosity
  grad_pz::Float64, # presure gradient
  Ha::Float64,      # Hartmann number
  n::Int,           # number of sumands included in Fourier series
  x)                # evaluation point

  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  V = 0.0; V0=0.0;
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = (Ha^2 + 4*α_k^2)^(0.5)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η))+exp(-r1_k*(1+η))
    den = 1+exp(-2*r1_k)
    V2 = (r2_k/N)*(num/den)

    num = exp(-r2_k*(1-η))+exp(-r2_k*(1+η))
    den = 1+exp(-2*r2_k)
    V3 = (r1_k/N)*(num/den)


    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3) * (1-V2-V3)
  end
  u_z = V/μ * (-grad_pz) * a^2

  VectorValue(0.0*u_z,0.0*u_z,u_z)
end


function analytical_hunt_j(
  a::Float64,       # semi-length of side walls
  b::Float64,       # semi-length of Hartmann walls
  σ::Float64,       # fluid conductivity
  μ::Float64,       # fluid viscosity
  grad_pz::Float64, # presure gradient
  Ha::Float64,      # Hartmann number
  n::Int,           # number of sumands included in Fourier series
  x)                # evaluation point

  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  H_dx = 0.0; H_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = sqrt(Ha^2 + 4*α_k^2)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η)) - exp(-r1_k*(1+η))
    num_dy = exp(-r1_k*(1-η))*(r1_k/a) + exp(-r1_k*(1+η))*(r1_k/a)
    den = 1+exp(-2*r1_k)
    H2 = (r2_k/N)*(num/den)
    H2_dy = (r2_k/N)*(num_dy/den)

    num = exp(-r2_k*(1-η)) - exp(-r2_k*(1+η))
    num_dy = exp(-r2_k*(1-η))*(r2_k/a) + exp(-r2_k*(1+η))*(r2_k/a)
    den = 1+exp(-2*r2_k)
    H3 = (r1_k/N)*(num/den)
    H3_dy = (r1_k/N)*(num_dy/den)

    H_dx += -2*(-1)^k * sin(α_k * ξ)/(a*l*α_k^2) * (H2 - H3)
    H_dy += 2*(-1)^k * cos(α_k * ξ)/(l*α_k^3) * (H2_dy - H3_dy)
  end
  j_x = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * H_dy
  j_y = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * (-H_dx)

  VectorValue(j_x,j_y,0.0)
end



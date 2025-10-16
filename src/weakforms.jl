
function weak_form(params)
  formulation = params[:fespaces][:formulation]
  if formulation == :H1HDiv
    if has_transient(params)
      weak_form_h1_hdiv_transient(params) 
    else
      weak_form_h1_hdiv(params)
    end
  elseif formulation == :H1H1
    if has_transient(params)
      weak_form_h1_h1_transient(params)
    else
      weak_form_h1_h1(params)
    end
  elseif formulation == :HDivHDiv
    if has_transient(params)
      weak_form_hdiv_hdiv_transient(params)
    else
      weak_form_hdiv_hdiv(params)
    end
  elseif formulation == :HDivH1
    if has_transient(params)
      weak_form_hdiv_h1_transient(params)
    else
      weak_form_hdiv_h1(params)
    end
  else
    @error("Unsupported formulation: $formulation")
  end
end

############################################################################################
# Variable management

setup_variable(x) = setup_variable(x...)

# H1-HDiv and HDiv-HDiv formulations
function setup_variable(u,p,j,φ)
  ∇u, ∇j = ∇(u), ∇(j)
  divu, divj = Operation(tr)(∇u), Operation(tr)(∇j)
  return (; u, p, j, φ, ∇u, divu, ∇j, divj)
end

# H1-H1 and HDiv-H1 formulation
function setup_variable(u,p,φ)
  ∇u, ∇φ = ∇(u), ∇(φ)
  divu = Operation(tr)(∇u)
  return (; u, p, φ, ∇u, divu, ∇φ)
end

# U-block GMG solver
function setup_variable_u(u)
  ∇u = ∇(u)
  divu = Operation(tr)(∇u)
  return (; u, ∇u, divu)
end

# UJ-block GMG solver
setup_variable_uj(x) = setup_variable_uj(x...)

function setup_variable_uj(u,j)
  ∇u, ∇j = ∇(u), ∇(j)
  divu, divj = Operation(tr)(∇u), Operation(tr)(∇j)
  return (; u, j, ∇u, divu, ∇j, divj)
end

############################################################################################
# Parameter retrieval

retrieve_fluid_params(params) = retrieve_fluid_params(params[:model],params)

function retrieve_fluid_params(model,params)
  fluid = params[:fluid]
  Ωf  = interior(params,model,fluid[:domain])
  dΩf = measure(params,Ωf)

  α, β, γ, σf = fluid[:α], fluid[:β], fluid[:γ], fluid[:σ]
  f, B, ζᵤ, ζⱼ = fluid[:f], fluid[:B], fluid[:ζᵤ], fluid[:ζⱼ]
  g, divg = fluid[:g], fluid[:divg]
  Πp = local_projection_operator(params)
  return α, β, γ, B, σf, f, g, divg, ζᵤ, ζⱼ, Πp, fluid[:convection], dΩf
end

retrieve_hdiv_fluid_params(params) = retrieve_hdiv_fluid_params(params[:model],params)

function retrieve_hdiv_fluid_params(model,params)
  Ωf = interior(params,model,params[:fluid][:domain])
  μ = params[:fluid][:μ]

  Γ = boundary(params,Ωf,nothing)
  h_Γ = get_cell_size(Γ)
  n_Γ = normal_vector(params,Γ)
  dΓ = measure(params,Γ)
  
  Λ = skeleton(params,Ωf,nothing)
  h_Λ = get_cell_size(Λ)
  n_Λ = normal_vector(params,Λ)
  dΛ = measure(params,Λ)

  if isa(params[:bcs][:u][:tags],Array)
    tags = params[:bcs][:u][:tags]
    values = params[:bcs][:u][:values]
  else
    tags = [params[:bcs][:u][:tags]]
    values = [params[:bcs][:u][:values]]
  end

  ΓD_params = []
  for (tag,u_D) in zip(tags,values)
    Γ_D = boundary(params,Ωf,tag)
    dΓ_D = measure(params,Γ_D)
    n_Γ_D = normal_vector(params,Γ_D)
    h_Γ_D = get_cell_size(Γ_D)
    push!(ΓD_params,(u_D,n_Γ_D,h_Γ_D,dΓ_D))
  end
  
  return μ,h_Γ,n_Γ,dΓ,h_Λ,n_Λ,dΛ,ΓD_params
end

retrieve_solid_params(params) = retrieve_solid_params(params[:model],params)

function retrieve_solid_params(model,params)
  solid  = params[:solid]
  if has_solid(params)
    Ωs  = interior(params,model,solid[:domain])
    dΩs = measure(params,Ωs)
    σs  = solid[:σ]
    ζ = solid[:ζ]
    g = solid[:g]
    divg = solid[:divg]
    return σs, g, divg, ζ, dΩs
  end
  return nothing
end

retrieve_bcs_params(params) = retrieve_bcs_params(params[:model],params)

function retrieve_bcs_params(model,params)
  bcs = params[:bcs]

  params_φ = []
  if isa(bcs[:φ],Array)
    for i in 1:length(bcs[:φ])
      φ_i = bcs[:φ][i][:value]
      Γ   = boundary(params,bcs[:φ][i][:domain])
      dΓ  = measure(params,Γ)
      n_Γ = normal_vector(params,Γ)
      push!(params_φ,(φ_i,n_Γ,dΓ))
    end
  end

  params_thin_wall = []
  for i in 1:length(bcs[:thin_wall])
    τ_i  = bcs[:thin_wall][i][:τ]
    cw_i = bcs[:thin_wall][i][:cw]
    jw_i = bcs[:thin_wall][i][:jw]
    Γ    = boundary(params,bcs[:thin_wall][i][:domain])
    dΓ   = measure(params,Γ)
    n_Γ  = normal_vector(params,Γ)
    push!(params_thin_wall,(τ_i,cw_i,jw_i,n_Γ,dΓ))
  end

  if length(bcs[:t]) != 0
    error("Boundary tranction not yet implemented")
  end

  params_Λ = []
  for i in 1:length(params[:bcs][:stabilization])
    Λ = skeleton(params,params[:bcs][:stabilization][i][:domain])
    dΛ = measure(params,Λ)
    h = get_cell_size(Λ)
    μ = params[:bcs][:stabilization][i][:μ]
    push!(params_Λ,(μ,h,dΛ))
  end

  return params_φ, params_thin_wall, params_Λ
end

############################################################################################
# H1-HDiv formulation

weak_form_h1_hdiv(params) = weak_form_h1_hdiv(params[:model],params)

function weak_form_h1_hdiv(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...
  )
  res(x,dy) = res_h1_hdiv(x,dy,weakform_params)
  jac(x,dx,dy) = jac_h1_hdiv(x,dx,dy,weakform_params)
  return res, jac
end

weak_form_h1_hdiv_transient(params) = weak_form_h1_hdiv_transient(params[:model],params)

function weak_form_h1_hdiv_transient(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...
  )
  dΩf = last(first(weakform_params))
  res(t,x,dy) = res_h1_hdiv(x,dy,time_eval(weakform_params,t)) + res_transient(x,dy,dΩf)
  jac(t,x,dx,dy) = jac_h1_hdiv(x,dx,dy,time_eval(weakform_params,t))
  jac_t(t,x,dx,dy) = jac_transient(dx,dy,dΩf)
  return res, jac, jac_t
end

function res_h1_hdiv(_x, _dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ = params

  x = setup_variable(_x)
  dy = setup_variable(_dy)

  r = res_fluid_h1_hdiv(x,dy,fluid_params...)
  if !isnothing(solid_params)
    r = r + res_solid_h1_hdiv(x,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + res_thin_wall(x,dy,p...)
  end
  for p in params_φ
    r = r + res_φ_bcs(x,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dy,p...)
  end

  return r
end

function jac_h1_hdiv(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ = params

  x = setup_variable(_x)
  dx = setup_variable(_dx)
  dy = setup_variable(_dy)

  r = jac_fluid_h1_hdiv(x,dx,dy,fluid_params...)
  if !isnothing(solid_params)
    r = r + jac_solid_h1_hdiv(x,dx,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + jac_thin_wall(x,dx,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

function res_fluid_h1_hdiv(x,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, v = x[:u], dy[:u]
  p, q = x[:p], dy[:p]
  j, s = x[:j], dy[:j]
  φ, ϕ = x[:φ], dy[:φ]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  div_u, div_v = x[:divu], dy[:divu]
  div_j, div_s = x[:divj], dy[:divj]

  u_block = β*(∇u⊙∇v)
  j_block = j⋅s

  # Augmented Lagrangian terms
  if !iszero(ζᵤ)
    u_block += ζᵤ*(Πp(u)*div_v)
  end
  if !iszero(ζⱼ)
    j_block += ζⱼ*(div_j*div_s)
  end

  # Convection term
  if convection != :none
    u_block += α*v⋅(conv∘(u,∇u))
  end

  return ∫(u_block - p*div_v - γ*(j×B)⋅v - div_u*q + j_block - σ*φ*div_s - σ*(u×B)⋅s - div_j*ϕ - f⋅v - g⋅s) * dΩ
end

function jac_fluid_h1_hdiv(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, ∇u = x[:u], x[:∇u]
  du, v = dx[:u], dy[:u]
  dp, q = dx[:p], dy[:p]
  dj, s = dx[:j], dy[:j]
  dφ, ϕ = dx[:φ], dy[:φ]
  ∇du, ∇v = dx[:∇u], dy[:∇u]
  div_du, div_v = dx[:divu], dy[:divu]
  div_dj, div_s = dx[:divj], dy[:divj]

  u_block = β*(∇du⊙∇v)
  j_block = dj⋅s 

  # Augmented Lagrangian terms
  if !iszero(ζᵤ)
    u_block += ζᵤ*(Πp(du)*div_v)
  end
  if !iszero(ζⱼ)
    j_block += ζⱼ*(div_dj*div_s)
  end

  # Convection term
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block - dp*div_v - γ*(dj×B)⋅v - div_du*q + j_block - σ*dφ*div_s - σ*(du×B)⋅s - div_dj*ϕ)dΩ
end

function res_solid_h1_hdiv(x,dy,σ,g,divg,ζ,dΩ)
  j, s = x[:j], dy[:j]
  φ, ϕ = x[:φ], dy[:φ]
  div_j, div_s = x[:divj], dy[:divj]

  j_block = j⋅s
  if !iszero(ζ)
    j_block += ζ*(div_j*div_s)
  end

  return ∫(j_block - σ*φ*div_s + ϕ*div_j - s⋅g)*dΩ
end

function jac_solid_h1_hdiv(x,dx,dy,σ,g,divg,ζ,dΩ)
  j, s = dx[:j], dy[:j]
  φ, ϕ = dx[:φ], dy[:φ]
  div_j, div_s = dx[:divj], dy[:divj]

  j_block = j⋅s
  if !iszero(ζ)
    j_block += ζ*(div_j*div_s)
  end

  return ∫(j_block - σ*φ*div_s + ϕ*div_j)*dΩ
end

############################################################################################
# H1-H1 formulation

weak_form_h1_h1(params) = weak_form_h1_h1(params[:model],params)

function weak_form_h1_h1(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...
  )
  res(x,dy) = res_h1_h1(x,dy,weakform_params)
  jac(x,dx,dy) = jac_h1_h1(x,dx,dy,weakform_params)
  return res, jac
end

weak_form_h1_h1_transient(params) = weak_form_h1_h1_transient(params[:model],params)

function weak_form_h1_h1_transient(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...
  )
  dΩf = last(first(weakform_params))
  res(t,x,dy) = res_h1_h1(x,dy,time_eval(weakform_params,t)) + res_transient(x,dy,dΩf)
  jac(t,x,dx,dy) = jac_h1_h1(x,dx,dy,time_eval(weakform_params,t))
  jac_t(t,x,dx,dy) = jac_transient(dx,dy,dΩf)
  return res, jac, jac_t
end

function res_h1_h1(_x, _dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ = params

  x = setup_variable(_x)
  dy = setup_variable(_dy)

  r = res_fluid_h1_h1(x,dy,fluid_params...)
  if !isnothing(solid_params)
    r = r + res_solid_h1_h1(x,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + res_thin_wall(x,dy,p...)
  end
  for p in params_φ
    r = r + res_φ_bcs(x,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dy,p...)
  end

  return r
end

function jac_h1_h1(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ = params

  x = setup_variable(_x)
  dx = setup_variable(_dx)
  dy = setup_variable(_dy)

  r = jac_fluid_h1_h1(x,dx,dy,fluid_params...)
  if !isnothing(solid_params)
    r = r + jac_solid_h1_h1(x,dx,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + jac_thin_wall(x,dx,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

function res_fluid_h1_h1(x,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, v = x[:u], dy[:u]
  p, q = x[:p], dy[:p]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  div_u, div_v = x[:divu], dy[:divu]
  φ,  ϕ = x[:φ], dy[:φ]
  ∇φ, ∇ϕ = x[:∇φ], dy[:∇φ]

  uB, vB = u×B, v×B

  u_block = β*(∇u⊙∇v) + γ*uB⋅vB

  # Augmented Lagrangian term
  if !iszero(ζᵤ)
    u_block += ζᵤ*(Πp(u)*div_v)
  end

  # Convection term
  if convection != :none
    u_block += α*v⋅(conv∘(u,∇u))
  end

  return ∫(u_block - p*div_v - div_u*q + ∇φ⋅∇ϕ - γ*(∇φ⋅vB) - uB⋅∇ϕ - f⋅v - divg*ϕ) * dΩ
end

function jac_fluid_h1_h1(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, ∇u = x[:u], x[:∇u]
  du, v = dx[:u], dy[:u]
  dp, q = dx[:p], dy[:p]
  dφ, ϕ = dx[:φ], dy[:φ]
  ∇du, ∇v = dx[:∇u], dy[:∇u]
  div_du, div_v = dx[:divu], dy[:divu]
  ∇dφ, ∇ϕ = dx[:∇φ], dy[:∇φ]

  duB, vB = du×B, v×B

  u_block = β*(∇du⊙∇v) + γ*duB⋅vB

  # Augmented Lagrangian term
  if !iszero(ζᵤ)
    u_block += ζᵤ*(Πp(du)*div_v)
  end

  # Convection term
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block - dp*div_v - div_du*q + ∇dφ⋅∇ϕ - γ*(∇dφ⋅vB) - duB⋅∇ϕ) * dΩ
end

function res_solid_h1_h1(x,dy,σ,g,divg,ζ,dΩ)
  φ, ϕ = x[:φ], dy[:φ]
  ∇φ, ∇ϕ = x[:∇φ], dy[:∇φ]
  return ∫(∇φ⋅∇ϕ - divg*ϕ)*dΩ
end

function jac_solid_h1_h1(x,dx,dy,σ,g,divg,ζ,dΩ)
  dφ, ϕ = dx[:φ], dy[:φ]
  ∇dφ, ∇ϕ = dx[:∇φ], dy[:∇φ]
  return ∫(∇dφ⋅∇ϕ)*dΩ
end

############################################################################################
# HDiv - HDiv formulation

weak_form_hdiv_hdiv(params) = weak_form_hdiv_hdiv(params[:model],params)

function weak_form_hdiv_hdiv(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)..., 
    retrieve_hdiv_fluid_params(model,params)
  )
  res(x,dy) = res_hdiv_hdiv(x,dy,weakform_params)
  jac(x,dx,dy) = jac_hdiv_hdiv(x,dx,dy,weakform_params)
  return res, jac
end

weak_form_hdiv_hdiv_transient(params) = weak_form_hdiv_hdiv_transient(params[:model],params)

function weak_form_hdiv_hdiv_transient(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)..., 
    retrieve_hdiv_fluid_params(model,params)
  )
  dΩf = last(first(weakform_params))
  res(t,x,dy) = res_hdiv_hdiv(x,dy,time_eval(weakform_params,t)) + res_transient(x,dy,dΩf)
  jac(t,x,dx,dy) = jac_hdiv_hdiv(x,dx,dy,time_eval(weakform_params,t))
  jac_t(t,x,dx,dy) = jac_transient(dx,dy,dΩf)
  return res, jac, jac_t
end

function res_hdiv_hdiv(_x, _dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ, hdiv_params = params

  x = setup_variable(_x)
  dy = setup_variable(_dy)

  r = res_fluid_h1_hdiv(x,dy,fluid_params...)
  r = r + res_fluid_hdiv_stab(x,dy,hdiv_params...)
  if !isnothing(solid_params)
    r = r + res_solid_h1_hdiv(x,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + res_thin_wall(x,dy,p...)
  end
  for p in params_φ
    r = r + res_φ_bcs(x,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dy,p...)
  end

  return r
end

function jac_hdiv_hdiv(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ, hdiv_params = params

  x = setup_variable(_x)
  dx = setup_variable(_dx)
  dy = setup_variable(_dy)

  r = jac_fluid_h1_hdiv(x,dx,dy,fluid_params...)
  r = r + jac_fluid_hdiv_stab(x,dx,dy,hdiv_params...)
  if !isnothing(solid_params)
    r = r + jac_solid_h1_hdiv(x,dx,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + jac_thin_wall(x,dx,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

function res_fluid_hdiv_stab(x,dy,μ,h_Γ,n_Γ,dΓ,h_Λ,n_Λ,dΛ,ΓD_params)
  u, v = x[:u], dy[:u]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  uᵗ, vᵗ = jump(u⊗n_Λ), jump(v⊗n_Λ)
  αΛ, αΓ = μ/h_Λ, μ/h_Γ

  c  = ∫(αΛ*uᵗ⊙vᵗ - vᵗ⊙mean(∇u) - mean(∇v)⊙uᵗ)dΛ
  c += ∫(αΓ*u⋅v - v⋅(n_Γ⋅∇u) - (n_Γ⋅∇v)⋅u)dΓ

  for (u_D, n_Γ_D, h_Γ_D, dΓ_D) in ΓD_params
    αΓD = μ/h_Γ_D
    c -= ∫(αΓD*v⋅u_D - (n_Γ_D⋅∇v)⋅u_D)*dΓ_D
  end

  return c
end

function jac_fluid_hdiv_stab(x,dx,dy,μ,h_Γ,n_Γ,dΓ,h_Λ,n_Λ,dΛ,ΓD_params)
  u, v = dx[:u], dy[:u]
  ∇u, ∇v = dx[:∇u], dy[:∇u]
  uᵗ, vᵗ = jump(u⊗n_Λ), jump(v⊗n_Λ)
  αΛ, αΓ = μ/h_Λ, μ/h_Γ

  c  = ∫( αΛ*uᵗ⊙vᵗ - vᵗ⊙mean(∇u) - mean(∇v)⊙uᵗ)dΛ
  c += ∫(αΓ*u⋅v - v⋅(n_Γ⋅∇u) - (n_Γ⋅∇v)⋅u)dΓ

  return c
end

############################################################################################
# HDiv-H1 formulation

weak_form_hdiv_h1(params) = weak_form_hdiv_h1(params[:model],params)

function weak_form_hdiv_h1(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...,
    retrieve_hdiv_fluid_params(model,params)
  )
  res(x,dy) = res_hdiv_h1(x,dy,weakform_params)
  jac(x,dx,dy) = jac_hdiv_h1(x,dx,dy,weakform_params)
  return res, jac
end

weak_form_hdiv_h1_transient(params) = weak_form_hdiv_h1_transient(params[:model],params)

function weak_form_hdiv_h1_transient(model,params)
  weakform_params = (
    retrieve_fluid_params(model,params), 
    retrieve_solid_params(model,params), 
    retrieve_bcs_params(model,params)...,
    retrieve_hdiv_fluid_params(model,params)
  )
  dΩf = last(first(weakform_params))
  res(t,x,dy) = res_hdiv_h1(x,dy,time_eval(weakform_params,t)) + res_transient(x,dy,dΩf)
  jac(t,x,dx,dy) = jac_hdiv_h1(x,dx,dy,time_eval(weakform_params,t))
  jac_t(t,x,dx,dy) = jac_transient(dx,dy,dΩf)
  return res, jac, jac_t
end

function res_hdiv_h1(_x, _dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ, hdiv_params = params

  x = setup_variable(_x)
  dy = setup_variable(_dy)

  r = res_fluid_h1_h1(x,dy,fluid_params...)
  r = r + res_fluid_hdiv_stab(x,dy,hdiv_params...)
  if !isnothing(solid_params)
    r = r + res_solid_h1_h1(x,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + res_thin_wall(x,dy,p...)
  end
  for p in params_φ
    r = r + res_φ_bcs(x,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dy,p...)
  end

  return r
end

function jac_hdiv_h1(_x,_dx,_dy, params)
  fluid_params, solid_params, params_φ, params_thin_wall, params_Λ, hdiv_params = params

  x = setup_variable(_x)
  dx = setup_variable(_dx)
  dy = setup_variable(_dy)

  r = jac_fluid_h1_h1(x,dx,dy,fluid_params...)
  r = r + jac_fluid_hdiv_stab(x,dx,dy,hdiv_params...)
  if !isnothing(solid_params)
    r = r + jac_solid_h1_h1(x,dx,dy,solid_params...)
  end
  for p in params_thin_wall
    r = r + jac_thin_wall(x,dx,dy,p...)
  end
  for p in params_Λ
    r = r + a_Λ(x,dx,p...)
  end

  return r
end

############################################################################################
# Utils 

conv(u,∇u) = (∇u')⋅u

function local_projection_operator(params)
  # If pressure-robust, no need to project
  fluid_is_pressure_robust(params) && return divergence
  
  # Otherwise: 
  reffe_p = params[:fespaces][:reffe_p]
  qdegree = params[:fespaces][:q]
  Πp = MultilevelTools.LocalProjectionMap(divergence,reffe_p,qdegree)
  return Πp
end

# Skeleton stabilisation

function a_Λ(x,y,μ,h,dΛ)
  ∇u, ∇v = x[:∇u], y[:∇u]
  ∫( (1/2) * μ * (h*h) * jump( ∇u ) ⊙ jump( ∇v ))*dΛ
end

# Augmented Lagrangian

function res_al_u(x,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, div_v = x[:u], dy[:divu]
  return ∫(ζᵤ*(Πp(u)*div_v)) * dΩ
end

function jac_al_u(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  du, div_dv = dx[:u], dy[:divu]
  return ∫(ζᵤ*(Πp(du)*div_dv)) * dΩ
end

function res_al_uj(x,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  u, div_v = x[:u], dy[:divu]
  div_j, div_s = x[:divj], dy[:divj]
  return ∫(ζᵤ*(Πp(u)*div_v) + ζⱼ*(div_j*div_s)) * dΩ
end

function jac_al_uj(x,dx,dy,α,β,γ,B,σ,f,g,divg,ζᵤ,ζⱼ,Πp,convection,dΩ)
  du, div_v = dx[:u], dy[:divu]
  div_dj, div_s = dx[:divj], dy[:divj]
  return ∫(ζᵤ*(Πp(du)*div_v) + ζⱼ*(div_dj*div_s)) * dΩ
end

# Boundary conditions

function res_φ_bcs(x,dy,φ0,n_Γ,dΓ)
  s = dy[:j]
  return ∫( (s⋅n_Γ)*φ0 )*dΓ
end

function res_thin_wall(x,dy,τ,cw,jw,n_Γ,dΓ)
  jn = dx[:j]⋅n_Γ # Normal component of j
  sn = dy[:j]⋅n_Γ # Normal component of s
  ∇jnn = n_Γ⋅(dx[:∇j]⋅n_Γ) # Normal-Normal component of ∇j
  return ∫(τ*sn*(jn - jw + cw*∇jnn))*dΓ
end

function jac_thin_wall(x,dx,dy,τ,cw,jw,n_Γ,dΓ)
  jn = dx[:j]⋅n_Γ # Normal component of j
  sn = dy[:j]⋅n_Γ # Normal component of s
  ∇jnn = n_Γ⋅(dx[:∇j]⋅n_Γ) # Normal-Normal component of ∇j
  return ∫(τ*sn*(jn + cw*∇jnn))*dΓ
end

# Transient

function res_transient(x,dy,dΩ)
  u, v = first(x), first(dy)
  return ∫(∂t(u)⋅v )*dΩ
end

function jac_transient(x,dy,dΩ)
  u, v = first(x), first(dy)
  return ∫(u⋅v)*dΩ
end

time_eval(a::AbstractVector,t::Real) = map(f->time_eval(f,t),a)
time_eval(a::Tuple,t::Real) = map(f->time_eval(f,t),a)
time_eval(a::Function,t::Real) = a(t)
time_eval(a::Any,t::Real) = a

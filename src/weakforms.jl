
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

# H1-H1 formulation
function setup_variable(u,p,φ)
  ∇u, ∇φ = ∇(u), ∇(φ)
  divu = Operation(tr)(∇u)
  return (; u, p, φ, ∇u, divu, ∇φ)
end

############################################################################################
# Parameter retrieval

retrieve_fluid_params(params) = retrieve_fluid_params(params[:model],params)

function retrieve_fluid_params(model,params)
  fluid = params[:fluid]
  Ωf  = params[:Ωf]
  dΩf = measure(params,Ωf)

  α, β, γ, σf = fluid[:α], fluid[:β], fluid[:γ], fluid[:σ]
  f, B, ζ, g = fluid[:f], fluid[:B], fluid[:ζ], fluid[:g]
  Πp = local_projection_operator(params)
  return α, β, γ, B, σf, f, g, ζ, Πp, fluid[:convection], dΩf
end

retrieve_hdiv_fluid_params(params) = retrieve_hdiv_fluid_params(params[:model],params)

function retrieve_hdiv_fluid_params(model,params)
  Ωf  = params[:Ωf]

  Γ = boundary(params,Ωf,nothing)
  Λ = skeleton(params,Ωf,nothing)
  Γ_D = boundary(params,Ωf,params[:bcs][:u][:tags])

  h_Γ = get_cell_size(Γ)
  h_Λ = get_cell_size(Λ)
  n_Γ_D = normal_vector(params,Γ_D)
  n_Λ = normal_vector(params,Λ)

  dΓ = measure(params,Γ)
  dΓ_D = measure(params,Γ_D)
  dΛ = measure(params,Λ)

  μ = 100.0
  u_D = params[:bcs][:u][:values]

  return μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ
end

retrieve_solid_params(params) = retrieve_solid_params(params[:model],params)

function retrieve_solid_params(model,params)
  solid  = params[:solid]
  if has_solid(params)
    Ωs  = params[:Ωs]
    dΩs = measure(params,Ωs)
    σs  = solid[:σ]
    return σs, dΩs
  end
  return nothing
end

retrieve_bcs_params(params) = retrieve_bcs_params(params[:model],params)

function retrieve_bcs_params(model,params)
  bcs = params[:bcs]

  params_φ = []
  for i in 1:length(bcs[:φ])
    φ_i = bcs[:φ][i][:value]
    Γ   = boundary(params,bcs[:φ][i][:domain])
    dΓ  = measure(params,Γ)
    n_Γ = normal_vector(params,Γ)
    push!(params_φ,(φ_i,n_Γ,dΓ))
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

function weak_form_h1_hdiv(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)...
  )
  res(x,dy) = res_h1_hdiv(x,dy,weakform_params)
  jac(x,dx,dy) = jac_h1_hdiv(x,dx,dy,weakform_params)
  return res, jac
end

function weak_form_h1_hdiv_transient(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)...
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

function res_fluid_h1_hdiv(x,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
  u, v = x[:u], dy[:u]
  p, q = x[:p], dy[:p]
  j, s = x[:j], dy[:j]
  φ, ϕ = x[:φ], dy[:φ]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  div_u, div_v = x[:divu], dy[:divu]
  div_j, div_s = x[:divj], dy[:divj]

  u_block = β*(∇u⊙∇v)
  j_block = j⋅s

  # Augmented Lagrangian term
  if !iszero(ζ)
    u_block += ζ*(Πp(div_u)*div_v)
    j_block += ζ*(div_j*div_s)
  end

  # Convection term
  if convection != :none
    u_block += α*v⋅(conv∘(u,∇u))
  end

  return ∫(u_block - p*div_v - γ*(j×B)⋅v - div_u*q + j_block - σ*φ*div_s - σ*(u×B)⋅s - div_j*ϕ - f⋅v - g⋅s) * dΩ
end

function jac_fluid_h1_hdiv(x,dx,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
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

  # Augmented Lagrangian term
  if !iszero(ζ)
    u_block += ζ*(Πp(div_du)*div_v)
    j_block += ζ*(div_dj*div_s)
  end

  # Convection term
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block - dp*div_v - γ*(dj×B)⋅v - div_du*q + j_block - σ*dφ*div_s - σ*(du×B)⋅s - div_dj*ϕ)dΩ
end

function res_solid_h1_hdiv(x,dy,σ,g,ζ,dΩ)
  j, s = x[:j], dy[:j]
  φ, ϕ = x[:φ], dy[:φ]
  div_j, div_s = x[:divj], dy[:divj]

  j_block = j⋅s
  if !iszero(ζ)
    j_block += ζ*(div_j*div_s)
  end

  return ∫(j_block - σ*φ*div_s + ϕ*div_j - s⋅g)*dΩ
end

function jac_solid_h1_hdiv(x,dx,dy,σ,g,ζ,dΩ)
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

function weak_form_h1_h1(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)...
  )
  res(x,dy) = res_h1_h1(x,dy,weakform_params)
  jac(x,dx,dy) = jac_h1_h1(x,dx,dy,weakform_params)
  return res, jac
end

function weak_form_h1_h1_transient(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)...
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

function res_fluid_h1_h1(x,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
  u, v = x[:u], dy[:u]
  p, q = x[:p], dy[:p]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  div_u, div_v = x[:divu], dy[:divu]
  ∇φ, ∇ϕ = x[:∇φ], dy[:∇φ]

  uB, vB = u×B, v×B

  u_block = β*(∇u⊙∇v) + γ*uB⋅vB

  # Augmented Lagrangian term
  if !iszero(ζ)
    u_block += ζ*(Πp(div_u)*div_v)
  end

  # Convection term
  if convection != :none
    u_block += α*v⋅(conv∘(u,∇u))
  end

  return ∫(u_block - p*div_v - div_u*q + ∇φ⋅∇ϕ - γ*(∇φ⋅vB) - uB⋅∇ϕ - f⋅v) * dΩ
end

function jac_fluid_h1_h1(x,dx,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
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
  if !iszero(ζ)
    u_block += ζ*(Πp(div_du)*div_v)
  end

  # Convection term
  if convection == :picard
    u_block += α*v⋅(conv∘(u,∇du))
  elseif convection == :newton
    u_block += α*v⋅(conv∘(u,∇du) + conv∘(du,∇u))
  end

  return ∫(u_block - dp*div_v - div_du*q + ∇dφ⋅∇ϕ - γ*(∇dφ⋅vB) - duB⋅∇ϕ) * dΩ
end

function res_solid_h1_h1(x,dy,σ,g,ζ,dΩ)
  φ, ϕ = x[:φ], dy[:φ]
  ∇φ, ∇ϕ = x[:∇φ], dy[:∇φ]
  return ∫(∇φ⋅∇ϕ)*dΩ
end

function jac_solid_h1_h1(x,dx,dy,σ,g,ζ,dΩ)
  dφ, ϕ = dx[:φ], dy[:φ]
  ∇dφ, ∇ϕ = dx[:∇φ], dy[:∇φ]
  return ∫(∇dφ⋅∇ϕ)*dΩ
end

############################################################################################
# HDiv - HDiv formulation

function weak_form_hdiv_hdiv(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)..., retrieve_hdiv_fluid_params(params)
  )
  res(x,dy) = res_hdiv_hdiv(x,dy,weakform_params)
  jac(x,dx,dy) = jac_hdiv_hdiv(x,dx,dy,weakform_params)
  return res, jac
end

function weak_form_hdiv_hdiv_transient(params)
  weakform_params = (
    retrieve_fluid_params(params), retrieve_solid_params(params), retrieve_bcs_params(params)..., retrieve_hdiv_fluid_params(params)
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
  r += res_fluid_hdiv_stab(x,dy,hdiv_params...)
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
  r += jac_fluid_hdiv_stab(x,dy,hdiv_params...)
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

function res_fluid_hdiv_stab(x,dy,μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ)
  u, v = x[:u], dy[:u]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  uᵗ, vᵗ = jump(u⊗n_Λ), jump(v⊗n_Λ)
  αΛ, αΓ = μ/h_Λ, μ/h_Γ

  c  = ∫( αΛ*vᵗ⊙uᵗ - vᵗ⊙mean(∇u) - mean(∇v)⊙uᵗ)dΛ
  c -= ∫(v⋅(∇u⋅n_Γ_D) + (∇v⋅n_Γ_D)⋅(u-u_D))dΓ_D
  c += ∫(αΓ*v⋅(u-u_D))dΓ

  return c
end

function jac_fluid_hdiv_stab(x,dy,μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ)
  u, v = x[:u], dy[:u]
  ∇u, ∇v = x[:∇u], dy[:∇u]
  uᵗ, vᵗ = jump(u⊗n_Λ), jump(v⊗n_Λ)
  αΛ, αΓ = μ/h_Λ, μ/h_Γ

  c  = ∫( αΛ*vᵗ⊙uᵗ - vᵗ⊙mean(∇u) - mean(∇v)⊙uᵗ)dΛ
  c -= ∫(v⋅(∇u⋅n_Γ_D) + (∇v⋅n_Γ_D)⋅u)dΓ_D
  c += ∫(αΓ*v⋅u)dΓ

  return c
end

############################################################################################
# Utils 

conv(u,∇u) = (∇u')⋅u

function local_projection_operator(params)
  poly = params[:fespaces][:poly]
  fluid_disc = params[:fespaces][:fluid_disc]

  # If pressure-robust, no need to project
  A = (poly == TET) && fluid_disc ∈ (:SV,:Pk_dPkm1)
  B = fluid_disc ∈ (:RT,:BDM)
  (A || B) && return identity
  
  # Otherwise: 
  reffe_p = params[:fespaces][:reffe_p]
  qdegree = params[:fespaces][:q]
  Πp = MultilevelTools.LocalProjectionMap(identity,reffe_p,qdegree)
  return Πp
end

# Skeleton stabilisation

function a_Λ(x,dy,μ,h,dΛ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( (1/2) * μ * (h*h) * jump( ∇(u) ) ⊙ jump( ∇(v_u) ))*dΛ
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

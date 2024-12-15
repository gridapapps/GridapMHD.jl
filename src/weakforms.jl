
function weak_form(params)
  k = params[:fespaces][:k]
  weak_form(params,k)
end

function weak_form(params,k)
  if !has_transient(params)
    _weak_form(params,k)
  else
    _ode_weak_form(params,k)
  end
end

function _weak_form(params,k)
  Ω = params[:Ω]
  dΩ = Measure(Ω,2*k)

  fluid = params[:fluid]
  dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params)
  dΩs, σs = retrieve_solid_params(params)
  bcs_params = retrieve_bcs_params(params)
  hdiv_params = retrieve_hdiv_fluid_params(params)
  params_φ, params_thin_wall, params_f, params_B, params_Λ = bcs_params

  Πp = local_projection_operator(params,k)

  function a(x,dy)
    r = a_mhd(x,dy,β,γ,B,σf,dΩf)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,p...)
    end
    for p in params_B
      r = r + a_B(x,dy,p...)
    end
    for p in params_Λ
      r = r + a_Λ(x,dy,p...)
    end
    if has_solid(params)
      r = r + a_solid(x,dy,σs,dΩs)
    end
    if !isnothing(hdiv_params)
      r = r + a_HDiv(x,dy,hdiv_params...)
    end
    if abs(ζ) > eps(typeof(ζ))
      r = r + a_al(x,dy,ζ,Πp,dΩf,dΩ)
    end
    r
  end

  function ℓ(dy)
    r = ℓ_mhd(dy,f,dΩf)
    for p in params_φ
      r = r + ℓ_φ(dy,p...)
    end
    for p in params_thin_wall
      r = r + ℓ_thin_wall(dy,p...)
    end
    for p in params_f
      r = r + ℓ_f(dy,p...)
    end
    if !isnothing(hdiv_params)
      r = r + ℓ_HDiv(dy,hdiv_params...)
    end
    r
  end

  function c(x,dy)
    r = c_mhd(x,dy,α,dΩf)
    r
  end

  function dc(x,dx,dy)
    if fluid[:convection] == :picard
      r = p_dc_mhd(x,dx,dy,α,dΩf)
    else
      r = n_dc_mhd(x,dx,dy,α,dΩf)
    end
    r
  end

  function res(x,dy)
    r = a(x,dy) - ℓ(dy)
    if has_convection(params)
      r = r + c(x,dy)
    end
    r
  end
  function jac(x,dx,dy)
    r = a(dx,dy)
    if has_convection(params)
      r = r + dc(x,dx,dy)
    end
    r
  end
  return res, jac
end


function _ode_weak_form(params,k)

  fluid = params[:fluid]
  dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params)
  dΩs, σs = retrieve_solid_params(params)
  bcs_params = retrieve_bcs_params(params)
  hdiv_params = retrieve_hdiv_fluid_params(params)
  params_φ, params_thin_wall, params_f, params_B, params_Λ = bcs_params

  Πp = local_projection_operator(params,k)

  m(t,x,dy) = m_u(x,dy,dΩf)

  a_dt(t,x,dy) = a_dut(x,dy,dΩf)

  function a(t,x,dy)
    r = a_mhd(x,dy,β,γ,time_eval(B,t),σf,dΩf)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,time_eval(p,t)...)
    end
    for p in params_B
      r = r + a_B(x,dy,time_eval(p,t)...)
    end
    for p in params_Λ
      r = r + a_Λ(x,dy,p...)
    end
    if has_solid(params)
      r = r + a_solid(x,dy,σs,dΩs)
    end
    if !isnothing(hdiv_params)
      r = r + a_HDiv(x,dy,hdiv_params...)
    end
    if abs(ζ) > eps(typeof(ζ))
      r = r + a_al(x,dy,ζ,Πp,dΩf,dΩ)
    end
    r
  end

  function ℓ(t,dy)
    r = ℓ_mhd(dy,time_eval(f,t),dΩf) + ℓ_fj(dy,time_eval(g,t),dΩf)
    for p in params_φ
      r = r + ℓ_φ(dy,time_eval(p,t)...)
    end
    for p in params_thin_wall
      r = r + ℓ_thin_wall(dy,time_eval(p,t)...)
    end
    for p in params_f
      r = r + ℓ_f(dy,time_eval(p,t)...)
    end
    if !isnothing(hdiv_params)
      r = r + ℓ_HDiv(x,dy,hdiv_params...)
    end
    r
  end

  function c(x,dy)
    r = c_mhd(x,dy,α,dΩf)
    r
  end

  function dc(x,dx,dy)
    if fluid[:convection] == :picard
      r = p_dc_mhd(x,dx,dy,α,dΩf)
    else
      r = n_dc_mhd(x,dx,dy,α,dΩf)
    end
    r
  end

  res(t,x,dy) = c(x,dy) + a_dt(t,x,dy) + a(t,x,dy) - ℓ(t,dy)
  jac(t,x,dx,dy) = dc(x,dx,dy) + a(t,dx,dy)
  jac_t(t,x,dx,dy) = m(t,dx,dy)

  return res,jac,jac_t
end

############################################################################################
# Parameter retrieval

retrieve_fluid_params(params) = retrieve_fluid_params(params[:model],params)

function retrieve_fluid_params(model,params)
  fluid  = params[:fluid]
  Ωf  = params[:Ωf]
  dΩf = measure(params,Ωf)

  α, β, γ, σf = fluid[:α], fluid[:β], fluid[:γ], fluid[:σ]
  f, B, ζ, g = fluid[:f], fluid[:B], fluid[:ζ], fluid[:g]
  return dΩf, α, β, γ, σf, f, B, ζ, g
end

retrieve_hdiv_fluid_params(params) = retrieve_hdiv_fluid_params(params[:model],params)

function retrieve_hdiv_fluid_params(model,params)
  Ωf  = params[:Ωf]

  if has_hdiv_fluid_disc(params)
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

    return (μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ)
  else
    return nothing
  end
  return hdiv_params
end

retrieve_solid_params(params) = retrieve_solid_params(params[:model],params)

function retrieve_solid_params(model,params)
  solid  = params[:solid]
  if solid !== nothing
    Ωs  = params[:Ωs]
    dΩs = measure(params,Ωs)
    σs  = solid[:σ]
    return dΩs, σs
  else
    return nothing, nothing
  end
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

  params_f = []
  for i in 1:length(bcs[:f])
    f_i  = bcs[:f][i][:value]
    Ω_i  = interior(params,bcs[:f][i][:domain])
    dΩ_i = measure(params,Ω_i)
    push!(params_f,(f_i,dΩ_i))
  end

  params_B = []
  for i in 1:length(bcs[:B])
    B_i  = bcs[:B][i][:value]
    Ω_i  = interior(params,bcs[:B][i][:domain])
    dΩ_i = measure(params,Ω_i)
    push!(params_f,(γ,B_i,dΩ_i))
  end

  params_Λ = []
  for i in 1:length(params[:bcs][:stabilization])
    Λ = skeleton(params,params[:bcs][:stabilization][i][:domain])
    dΛ = measure(params,Λ)
    h = get_cell_size(Λ)
    μ = params[:bcs][:stabilization][i][:μ]
    push!(params_Λ,(μ,h,dΛ))
  end

  return params_φ, params_thin_wall, params_f, params_B, params_Λ
end

############################################################################################
# Weakform blocks

# MHD equations

conv(u,∇u) = (∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du) + conv(du,∇u)

function a_mhd(x,dy,β,γ,B,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫(
    β*(∇(u)⊙∇(v_u)) - p*(∇⋅v_u) -(γ*(j×B)⋅v_u) -
    (∇⋅u)*v_p +
    j⋅v_j - σ*φ*(∇⋅v_j) - σ*(u×B)⋅v_j +
    - (∇⋅j)*v_φ ) * dΩ
end

a_mhd_u_u(u,v_u,β,dΩ)   = ∫( β*(∇(u)⊙∇(v_u)) )*dΩ
a_mhd_u_p(p,v_u,dΩ)     = ∫( -p*(∇⋅v_u) )*dΩ
a_mhd_u_j(j,v_u,γ,B,dΩ) = ∫( -γ*(j×B)⋅v_u )*dΩ
a_mhd_p_u(u,v_p,dΩ)     = ∫( -(∇⋅u)*v_p )*dΩ

a_mhd_j_u(u,v_j,σ,B,dΩ) = ∫( -σ*(u×B)⋅v_j )*dΩ
a_mhd_j_j(j,v_j,dΩ)     = ∫( j⋅v_j )*dΩ
a_mhd_j_φ(φ,v_j,σ,dΩ)   = ∫( -σ*φ*(∇⋅v_j) )*dΩ
a_mhd_φ_j(j,v_φ,dΩ)     = ∫( -(∇⋅j)*v_φ )*dΩ

function ℓ_mhd(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ℓ_mhd_u(v_u,f,dΩ)
end
ℓ_mhd_u(v_u,f,dΩ) = ∫( v_u⋅f )*dΩ

# Convection

function c_mhd(x,dy,α,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  c_mhd_u_u(u,v_u,α,dΩ)
end
c_mhd_u_u(u,v_u,α,dΩ) = ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ

dc_mhd(x,dx,dy,α,dΩ) = n_dc_mhd(x,dx,dy,α,dΩ) # Default to full Newton iteration

# Convection derivative: Newton iteration
function n_dc_mhd(x,dx,dy,α,dΩ)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  n_dc_mhd_u_u(u,du,v_u,α,dΩ)
end
n_dc_mhd_u_u(u,du,v_u,α,dΩ) = ∫( α*v_u⋅( dconv∘(du,∇(du),u,∇(u)) ) ) * dΩ

# Convection derivative: Picard iteration
function p_dc_mhd(x,dx,dy,α,dΩ)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  p_dc_mhd_u_u(u,du,v_u,α,dΩ)
end
p_dc_mhd_u_u(u,du,v_u,α,dΩ) = ∫( α*v_u⋅( conv∘(u,∇(du)) ) ) * dΩ

# Stabilisation

function a_Λ(x,dy,μ,h,dΛ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( (1/2) * μ * (h*h) * jump( ∇(u) ) ⊙ jump( ∇(v_u) ))*dΛ
end

# Augmented lagrangian

function a_al(x,dy,ζ,Πp,dΩf,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  if isnothing(Πp)
    a_al_sf(u,v_u,ζ,dΩf) + a_al_sf(j,v_j,ζ,dΩ)
  else
    a_al_sf(u,v_u,ζ,Πp,dΩf) + a_al_sf(j,v_j,ζ,dΩ)
  end
end
a_al_sf(x,y,ζ,Π,dΩ) = ∫(ζ*Π(x)*(∇⋅y))*dΩ
a_al_sf(x,y,ζ,dΩ) = ∫(ζ*(∇⋅x)*(∇⋅y))*dΩ

function local_projection_operator(params,k)
  poly = params[:fespaces][:poly]
  fluid_disc = params[:fespaces][:fluid_disc]

  # If pressure-robust, no need to project
  A = (poly == TET) && fluid_disc ∈ (:SV,:Pk_dPkm1)
  B = fluid_disc ∈ (:RT,:BDM)
  if A || B
    return nothing
  end
  
  # Otherwise: 
  reffe_p = params[:fespaces][:reffe_p]
  Πp = MultilevelTools.LocalProjectionMap(divergence,reffe_p,2*k)
  return Πp
end

# Div-Conforming laplacian terms (parts integration + tangent component penalty)

function a_HDiv(x,y,μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ)
  u, _, _, _ = x
  v, _, _, _ = y

  αΛ, αΓ = μ/h_Λ, μ/h_Γ
  ∇u, ∇v = ∇(u), ∇(v)
  uᵗ, vᵗ = jump(u⊗n_Λ), jump(v⊗n_Λ)

  c  = ∫( αΛ*vᵗ⊙uᵗ - vᵗ⊙mean(∇u) - mean(∇v)⊙uᵗ)dΛ
  c -= ∫(v⋅(∇u⋅n_Γ_D) + (∇v⋅n_Γ_D)⋅u)dΓ_D
  c += ∫(αΓ*v⋅u)dΓ
  return c
end

function ℓ_HDiv(dy,μ,h_Γ,h_Λ,n_Γ_D,n_Λ,u_D,dΓ,dΓ_D,dΛ)
  v, _, _, _ = dy
  αΓ = μ/h_Γ
  c  = ∫(αΓ*v⋅u_D)dΓ - ∫((∇(v)⋅n_Γ_D)⋅u_D)dΓ_D
  return c
end

# Solid equations

function a_solid(x,dy,σ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( j⋅v_j - σ*φ*(∇⋅v_j) + (∇⋅j)*v_φ)dΩ
end
a_solid_j_j(j,v_j,dΩ)   = ∫( j⋅v_j )*dΩ
a_solid_j_φ(φ,v_j,σ,dΩ) = ∫( -σ*φ*(∇⋅v_j) )*dΩ
a_solid_φ_j(j,v_φ,dΩ)   = ∫( (∇⋅j)*v_φ )*dΩ

# Boundary conditions

function ℓ_φ(dy,φ,n_Γ,dΓ)
  v_u, v_p, v_j, v_φ = dy
  ∫( -(v_j⋅n_Γ)*φ )*dΓ
end
ℓ_φ_j(φ,v_j,n_Γ,dΓ) = ∫( -(v_j⋅n_Γ)*φ )*dΓ

function ℓ_f(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_u⋅f )*dΩ
end
ℓ_f_u(f,v_u,dΩ) = ∫( v_u⋅f )*dΩ

function ℓ_fj(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_j⋅f )*dΩ
end

function ℓ_thin_wall(dy,τ,cw,jw,n_Γ,dΓ)
  v_u, v_p, v_j, v_φ = dy
  ℓ_thin_wall_j(v_j,τ,cw,jw,n_Γ,dΓ)
end
ℓ_thin_wall_j(v_j,τ,cw,jw,n_Γ,dΓ) = ∫( τ*(v_j⋅n_Γ)*jw ) * dΓ

function a_thin_wall(x,dy,τ,cw,jw,n_Γ,dΓ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  a_thin_wall_j_j(j,v_j,τ,cw,jw,n_Γ,dΓ)
end
a_thin_wall_j_j(j,v_j,τ,cw,jw,n_Γ,dΓ) = ∫( τ*((v_j⋅n_Γ)*(j⋅n_Γ) + cw*(v_j⋅n_Γ)*(n_Γ⋅(∇(j)⋅n_Γ))) )*dΓ

function a_B(x,dy,γ,B,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( -(γ*(j×B)⋅v_u) - (u×B)⋅v_j )*dΩ
end

# Mass matrix

function a_dut(x,dy,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫(∂t(u)⋅v_u )*dΩ
end

function m_u(x,dy,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( u⋅v_u )*dΩ
end

############################################################################################
# Helper functions

time_eval(a::AbstractVector,t::Real) = map(f->time_eval(f,t),a)
time_eval(a::Tuple,t::Real) = map(f->time_eval(f,t),a)
time_eval(a::Function,t::Real) = a(t)
time_eval(a::Any,t::Real) = a

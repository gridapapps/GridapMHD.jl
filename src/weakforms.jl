
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
  Ωf, dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params,k)

  solid = params[:solid]
  Ωs, dΩs, σs = retrieve_solid_params(params,k)

  bcs_params = retrieve_bcs_params(params,k)
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
  Ωf, dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params,k)

  solid = params[:solid]
  Ωs, dΩs, σs = retrieve_solid_params(params,k)

  bcs_params = retrieve_bcs_params(params,k)
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

retrieve_fluid_params(params,k) = retrieve_fluid_params(params[:model],params,k)

function retrieve_fluid_params(model,params,k)
  fluid  = params[:fluid]
  quad3D = params[:fespaces][:quadratures][3]
  Ωf  = params[:Ωf]
  dΩf = Measure(Ωf,quad3D)

  α, β, γ, σf = fluid[:α], fluid[:β], fluid[:γ], fluid[:σ]
  f = fluid[:f]
  B = fluid[:B]
  ζ = fluid[:ζ]
  g = fluid[:g]
  return Ωf, dΩf, α, β, γ, σf, f, B, ζ, g
end

retrieve_solid_params(params,k) = retrieve_solid_params(params[:model],params,k)

function retrieve_solid_params(model,params,k)
  quad3D = params[:fespaces][:quadratures][3]
  solid  = params[:solid]
  if solid !== nothing
    Ωs  = params[:Ωs]
    dΩs = Measure(Ωs,quad3D)
    σs  = solid[:σ]
    return Ωs, dΩs, σs
  else
    return nothing, nothing, nothing
  end
end

retrieve_bcs_params(params,k) = retrieve_bcs_params(params[:model],params,k)

function retrieve_bcs_params(model,params,k)
  quad3D = params[:fespaces][:quadratures][3]
  quad2D = params[:fespaces][:quadratures][2]
  bcs = params[:bcs]

  params_φ = []
  for i in 1:length(bcs[:φ])
    φ_i = bcs[:φ][i][:value]
    Γ   = _boundary(model,bcs[:φ][i][:domain])
    dΓ  = Measure(Γ,quad2D)
    n_Γ = get_normal_vector(Γ)
    push!(params_φ,(φ_i,n_Γ,dΓ))
  end

  params_thin_wall = []
  for i in 1:length(bcs[:thin_wall])
    τ_i  = bcs[:thin_wall][i][:τ]
    cw_i = bcs[:thin_wall][i][:cw]
    jw_i = bcs[:thin_wall][i][:jw]
    Γ    = _boundary(model,bcs[:thin_wall][i][:domain])
    dΓ   = Measure(Γ,quad2D)
    n_Γ  = get_normal_vector(Γ)
    push!(params_thin_wall,(τ_i,cw_i,jw_i,n_Γ,dΓ))
  end

  if length(bcs[:t]) != 0
    error("Boundary tranction not yet implemented")
  end

  params_f = []
  for i in 1:length(bcs[:f])
    f_i  = bcs[:f][i][:value]
    Ω_i  = _interior(model,bcs[:f][i][:domain])
    dΩ_i = Measure(Ω_i,quad3D)
    push!(params_f,(f_i,dΩ_i))
  end

  params_B = []
  for i in 1:length(bcs[:B])
    B_i  = bcs[:B][i][:value]
    Ω_i  = _interior(model,bcs[:B][i][:domain])
    dΩ_i = Measure(Ω_i,quad3D)
    push!(params_f,(γ,B_i,dΩ_i))
  end

  params_Λ = []
  for i in 1:length(params[:bcs][:stabilization])
    Λ = _skeleton(model,params[:bcs][:stabilization][i][:domain])
    dΛ = Measure(Λ,quad2D)
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

function a_Λ(x,dy,μ,h,dΛ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( (1/2) * μ * (h*h) * jump( ∇(u) ) ⊙ jump( ∇(v_u) ))*dΛ
end

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
  fluid_disc = params[:fespaces][:fluid_disc]
  if fluid_disc ∈ (:SV,)
    return nothing
  end
  
  reffe_p = params[:fespaces][:reffe_p]
  Πp = MultilevelTools.LocalProjectionMap(divergence,reffe_p,2*k)
  return Πp
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

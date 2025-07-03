
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
  
  fluid_params = retrieve_fluid_params(params)
  solid_params = retrieve_solid_params(params)
  params_φ, params_thin_wall, params_Λ  = retrieve_bcs_params(params)

  function res(_x,_dy)
    x = setup_variable(_x)
    dy = setup_variable(_dy)

    r = res_fluid(x,dy,fluid_params...)
    if has_solid(params)
      r = r + res_solid(x,dy,solid_params...)
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

  function jac(_x,_dx,_dy)
    x = setup_variable(_x)
    dx = setup_variable(_dx)
    dy = setup_variable(_dy)

    r = jac_fluid(x,dx,dy,fluid_params...)
    if has_solid(params)
      r = r + jac_solid(x,dx,dy,solid_params...)
    end
    for p in params_thin_wall
      r = r + jac_thin_wall(x,dx,dy,p...)
    end
    for p in params_Λ
      r = r + a_Λ(x,dx,p...)
    end
    return r
  end

  return res, jac
end


function _ode_weak_form(params,k)

  fluid = params[:fluid]
  dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params)
  dΩs, σs = retrieve_solid_params(params)
  bcs_params = retrieve_bcs_params(params)
  hdiv_params = retrieve_hdiv_fluid_params(params)
  params_φ, params_thin_wall, params_Λ = bcs_params

  Πp = local_projection_operator(params,k)

  m(t,x,dy) = m_u(x,dy,dΩf)

  a_dt(t,x,dy) = a_dut(x,dy,dΩf)

  function a(t,x,dy)
    r = a_mhd(x,dy,β,γ,time_eval(B,t),σf,dΩf)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,time_eval(p,t)...)
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
# Variable management

setup_variable(x) = setup_variable(x...)

function setup_variable(u,p,j,φ)
  ∇u, ∇j= ∇(u), ∇(j)
  divu, divj = Operation(tr)(∇u), Operation(tr)(∇j)
  return (; u, p, j, φ, ∇u, divu, ∇j, divj)
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
# Weakform blocks

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

# Fluid equations

function res_fluid(x,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
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

function jac_fluid(x,dx,dy,α,β,γ,B,σ,f,g,ζ,Πp,convection,dΩ)
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

# Skeleton stabilisation

function a_Λ(x,dy,μ,h,dΛ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( (1/2) * μ * (h*h) * jump( ∇(u) ) ⊙ jump( ∇(v_u) ))*dΛ
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

function res_solid(x,dy,σ,g,ζ,dΩ)
  j, s = x[:j], dy[:j]
  φ, ϕ = x[:φ], dy[:φ]
  div_j, div_s = x[:divj], dy[:divj]

  j_block = j⋅s
  if !iszero(ζ)
    j_block += ζ*(div_j*div_s)
  end

  return ∫(j_block - σ*φ*div_s + ϕ*div_j - s⋅g)*dΩ
end

function jac_solid(x,dx,dy,σ,g,ζ,dΩ)
  j, s = dx[:j], dy[:j]
  φ, ϕ = dx[:φ], dy[:φ]
  div_j, div_s = dx[:divj], dy[:divj]

  j_block = j⋅s
  if !iszero(ζ)
    j_block += ζ*(div_j*div_s)
  end

  return ∫(j_block - σ*φ*div_s + ϕ*div_j)*dΩ
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

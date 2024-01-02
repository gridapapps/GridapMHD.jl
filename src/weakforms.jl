
function weak_form(params,k)

  fluid = params[:fluid]
  Ωf, dΩf, α, β, γ, σf, f, B, ζ = retrieve_fluid_params(params,k)

  solid = params[:solid]
  Ωs, dΩs, σs = retrieve_solid_params(params,k)

  params_φ, params_thin_wall, params_f, params_B = retrieve_bcs_params(params,k)

  function a(x,dy)
    r = a_mhd(x,dy,β,γ,B,σf,dΩf)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,p...)
    end
    for p in params_B
      r = r + a_B(x,dy,p...)
    end
    if solid !== nothing
      r = r + a_solid(x,dy,σs,dΩs)
    end
    if abs(ζ) > eps(typeof(ζ))
      r = r + a_augmented_lagragian(x,dy,ζ,dΩf)
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
    r = dc_mhd(x,dx,dy,α,dΩf)
    r
  end

  res(x,dy) = c(x,dy) + a(x,dy) - ℓ(dy)
  jac(x,dx,dy) = dc(x,dx,dy) + a(dx,dy)

  return res, jac
end

############################################################################################
# Parameter retrieval 

function retrieve_fluid_params(params,k)
  fluid = params[:fluid]
  Ωf    = _interior(params[:model],fluid[:domain])
  dΩf   = Measure(Ωf,2*k)

  α, β, γ, σf = fluid[:α], fluid[:β], fluid[:γ], fluid[:σ]
  f = fluid[:f]
  B = fluid[:B]
  ζ = fluid[:ζ]
  return Ωf, dΩf, α, β, γ, σf, f, B, ζ
end

function retrieve_solid_params(params,k)
  solid = params[:solid]
  if solid !== nothing
    Ωs  = _interior(params[:model],solid[:domain])
    dΩs = Measure(Ωs,2*k)
    σs  = solid[:σ]
    return Ωs, dΩs, σs
  else
    return nothing, nothing, nothing
  end
end

function retrieve_bcs_params(params,k)
  bcs = params[:bcs]

  params_φ = []
  for i in 1:length(bcs[:φ])
    φ_i = bcs[:φ][i][:value]
    Γ   = _boundary(params[:model],bcs[:φ][i][:domain])
    dΓ  = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_φ,(φ_i,n_Γ,dΓ))
  end

  params_thin_wall = []
  for i in 1:length(bcs[:thin_wall])
    τ_i  = bcs[:thin_wall][i][:τ]
    cw_i = bcs[:thin_wall][i][:cw]
    jw_i = bcs[:thin_wall][i][:jw]
    Γ    = _boundary(params[:model],bcs[:thin_wall][i][:domain])
    dΓ   = Measure(Γ,2*k)
    n_Γ  = get_normal_vector(Γ)
    push!(params_thin_wall,(τ_i,cw_i,jw_i,n_Γ,dΓ))
  end

  if length(bcs[:t]) != 0
    error("Boundary tranction not yet implemented")
  end

  params_f = []
  for i in 1:length(bcs[:f])
    f_i  = bcs[:f][i][:value]
    Ω_i  = _interior(params[:model],bcs[:f][i][:domain])
    dΩ_i = Measure(Ω_i,2*k)
    push!(params_f,(f_i,dΩ_i))
  end

  params_B = []
  for i in 1:length(bcs[:B])
    B_i  = bcs[:B][i][:value]
    Ω_i  = _interior(params[:model],bcs[:B][i][:domain])
    dΩ_i = Measure(Ω_i,2*k)
    push!(params_f,(γ,B_i,dΩ_i))
  end

  return params_φ, params_thin_wall, params_f, params_B
end

############################################################################################
# Weakform blocks

# MHD equations

conv(u,∇u) = (∇u')⋅u

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

function c_mhd(x,dy,α,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  c_mhd_u_u(u,v_u,α,dΩ)
end
c_mhd_u_u(u,v_u,α,dΩ) = ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ

function dc_mhd(x,dx,dy,α,dΩ)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  dc_mhd_u_u(u,du,v_u,α,dΩ)
end
dc_mhd_u_u(u,du,v_u,α,dΩ) = ∫( α*v_u⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u))) ) ) * dΩ


# Augmented lagrangian

function a_al(x,dy,ζ,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  a_al_u_u(u,v_u,ζ,dΩ) + a_al_j_j(j,v_j,ζ,dΩ)
end
a_al_u_u(u,v_u,ζ,dΩ) = ∫( ζ*(∇⋅u)*(∇⋅v_u) ) * dΩ
a_al_j_j(j,v_j,ζ,dΩ) = ∫( ζ*(∇⋅j)*(∇⋅v_j) ) * dΩ

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

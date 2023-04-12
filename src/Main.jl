
# MHD equations
# ∇⋅u = 0
# u⋅∇(u) -ν*Δ(u) + (1/ρ)*∇(p) - (1/ρ)*(j×B) = (1/ρ)*f
# ∇⋅j = 0
# j + σ*∇(φ) - σ*(u×B) = 0
# solving for u,p,j,φ for a given B,ν,ρ,σ
#
# One can provide characteristic quantities
# u0,B0,L
#
# To introduce these change of variables
# u = u0*ū, B = B0*B̄, j = σ*u0*B0*j̄, φ = u0*B0*L*φ̄ (I am not sure about L in φ)
# p = ρ*u0^2*p̄, f = (ρ*u0^2/L)*f̄ (option 1)
# p = σ*u0*B0*L*p̄, f = (σ*u0*B0^2)*f̄ (option 2)
#
# To solve the following equations in a scaled domain according to L
# ∇⋅ū = 0
# ∇⋅j̄ = 0
# j̄ + ∇(φ̄) - ū×B̄ = 0
# ū⋅∇(ū) -(1/Re)*Δ(ū) + ∇(p̄) - N*(j̄×B̄) = f̄ (Option 1, CFD)
# (1/N)*ū⋅∇(ū) - (1/Ha^2)*Δ(ū) + ∇(p̄) - (j̄×B̄) = f̄ (Option 2,MHD)
#
# In order to account for both options, the code solves these equations
#
# ∇⋅ū = 0
# ∇⋅j̄ = 0
# j̄ + ∇(φ̄) - ū×B̄ = 0
# α*ū⋅∇(ū) - β*Δ(ū) + ∇(p̄) - γ*(j̄×B̄) = f̄
#
# α = 1, β = (1/Re), γ = N (option 1,CFD)
# α = (1/N), β = (1/Ha^2), γ = 1 (option 2,MHD)
#
# Some boundary conditions
# (for simplicity we drop the bars, i.e. ū is simply u)
# 
# Velocity bc
# u = u (imposed strongly)
#
# Charge bc
# j = j (imposed strongly in the normal direction to the boundary)
#
# Traction bC
# # n⋅∇(u) - p*n = t
#
# Potential bc
# φ = φ (Imposed weakly)
#
# Thin wall bc
# j⋅n + cw*n⋅∇(j)⋅n = jw (imposed weakly via a penalty of value τ)


function main(params::Dict)

  defaults = Dict(
    :solver=>NLSolver(show_trace=true,method=:newton),
    :solver_postpro=> (x->nothing),
    :matrix_type=>SparseMatrixCSC{Float64,Int},
    :vector_type=>Vector{Float64},
    :ptimer=>PTimer(get_part_ids(sequential,1)),
    :fluid=>Dict(
      :k=>2,
      :φ=>[],
      :t=>[],
      :thin_wall=>[]
    )
  )
  add_defaults!(params,defaults)

  t = params[:ptimer]
  tic!(t;barrier=true)
  k = params[:fluid][:k]
  Ω = params[:fluid][:domain]
  T = Float64
  D = num_cell_dims(Ω)
  reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
  reffe_p = ReferenceFE(lagrangian,T,k-1;space=:P)
  reffe_j = ReferenceFE(raviart_thomas,T,k-1)
  reffe_φ = ReferenceFE(lagrangian,T,k-1)

  # Test spaces
  V_u = TestFESpace(Ω,reffe_u;dirichlet_tags=params[:fluid][:u][:tags])
  V_p = TestFESpace(Ω,reffe_p)
  V_j = TestFESpace(Ω,reffe_j;dirichlet_tags=params[:fluid][:j][:tags])
  V_φ = TestFESpace(Ω,reffe_φ;conformity=:L2)
  V = MultiFieldFESpace([V_u,V_p,V_j,V_φ])

  # Trial spaces
  # TODO improve for parallel computations
  tic!(t;barrier=true)
  z = zero(VectorValue{D,Float64})
  u_bc = params[:fluid][:u][:values]
  j_bc = params[:fluid][:j][:values]
  U_u = u_bc == z ? V_u : TrialFESpace(V_u,u_bc)
  U_j = j_bc == z ? V_j : TrialFESpace(V_j,j_bc)
  U_p = TrialFESpace(V_p)
  U_φ = TrialFESpace(V_φ)
  U = MultiFieldFESpace([U_u,U_p,U_j,U_φ])
  toc!(t,"fe_spaces")

  tic!(t;barrier=true)
  if params[:debug]
    Random.seed!(1234)
    vt = get_vector_type(U)
    free_ids = get_free_dof_ids(U)
    free_vals = _rand(vt,free_ids)
    xh = FEFunction(U,free_vals)
    toc!(t,"solve")
  else
    res, jac = weak_form(params,k)
    Tm = params[:matrix_type]
    Tv = params[:vector_type]
    assem = SparseMatrixAssembler(Tm,Tv,U,V)
    op = FEOperator(res,jac,U,V,assem)
    xh = zero(U)
    if params[:solve]
       solver = params[:solver]
       xh,cache = solve!(xh,solver,op)
       solver_postpro = params[:solver_postpro]
       solver_postpro(cache)
       toc!(t,"solve")
    end
    if params[:only_assemble]
      tic!(t;barrier=true)
      r = residual(op,xh)
      toc!(t,"residual")
      tic!(t;barrier=true)
      j = jacobian(op,xh)
      toc!(t,"jacobian")
    end
  end
  xh
end

function _rand(vt::Type{<:Vector{T}},r::AbstractUnitRange) where T
  rand(T,length(r))
end

function _rand(vt::Type{<:PVector{T,A}},ids::PRange) where {T,A}
  values = map_parts(ids.partition) do partition
    Tv = eltype(A)
    _rand(Tv,1:num_lids(partition))
  end
  PVector(values,ids)
end

function add_defaults!(params,defaults)
  for (key,val) in defaults
    if !haskey(params,key)
      params[key] = val
    elseif isa(val,AbstractDict)
      @assert isa(params[key],AbstractDict)
      add_defaults!(params[key],val)
    end
  end
end

function weak_form(params,k)

  fluid = params[:fluid]

  Ω = Triangulation(fluid[:domain])
  dΩ = Measure(Ω,2*k)
  D = num_cell_dims(Ω)
  z = zero(VectorValue{D,Float64})

  α = fluid[:α]
  β = fluid[:β]
  γ = fluid[:γ]
  f = isa(fluid[:f],VectorValue) ? fluid[:f] : z
  B = isa(fluid[:B],VectorValue) ? fluid[:B] : z

  params_φ = []
  for i in 1:length(fluid[:φ])
    φ_i = fluid[:φ][i][:value]
    Γ = fluid[:φ][i][:domain]
    dΓ = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_φ,(φ_i,n_Γ,dΓ))
  end

  params_thin_wall = []
  for i in 1:length(fluid[:thin_wall])
    cw_i = fluid[:thin_wall][i][:cw]
    jw_i = fluid[:thin_wall][i][:jw]
    Γ = fluid[:thin_wall][i][:domain]
    dΓ = Measure(Γ,2*k)
    n_Γ = get_normal_vector(Γ)
    push!(params_thin_wall,(cw_i,jw_i,n_Γ,dΓ))
  end

  params_f = []
  if !isa(fluid[:f],VectorValue)
    for i in 1:length(fluid[:f])
      f_i = fluid[:f][i][:value]
      Ω_i = fluid[:f][i][:domain]
      dΩ_i = Measure(Ω_i,2*k)
      push!(params_f,(f_i,dΩ_i))
    end
  end

  params_B = []
  if !isa(fluid[:B],VectorValue)
    for i in 1:length(fluid[:B])
      B_i = fluid[:B][i][:value]
      Ω_i = fluid[:B][i][:domain]
      dΩ_i = Measure(Ω_i,2*k)
      push!(params_f,(γ,B_i,dΩ_i))
    end
  end

  function a(x,dy)
    r = a_mhd(x,dy,β,γ,B,dΩ)
    for p in params_thin_wall
      r = r + a_thin_wall(x,dy,p...)
    end
    for p in params_B
      r = r + a_B(x,dy,p...)
    end
    r
  end

  function ℓ(dy)
    r = ℓ_mhd(dy,f,dΩ)
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
    r = c_mhd(x,dy,α,dΩ)
    r
  end

  function dc(x,dx,dy)
    r = dc_mhd(x,dx,dy,α,dΩ)
    r
  end

  res(x,dy) = c(x,dy) + a(x,dy) - ℓ(dy)
  jac(x,dx,dy) = dc(x,dx,dy) + a(dx,dy)

  res, jac
end

conv(u,∇u) = (∇u')⋅u

function a_mhd(x,dy,β,γ,B,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫(
    β*(∇(u)⊙∇(v_u)) - p*(∇⋅v_u) -(γ*(j×B)⋅v_u) +
    (∇⋅u)*v_p +
    j⋅v_j - φ*(∇⋅v_j) - (u×B)⋅v_j +
    (∇⋅j)*v_φ ) * dΩ
end

function ℓ_mhd(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_u⋅f )*dΩ
end

function c_mhd(x,dy,α,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( α*v_u⋅(conv∘(u,∇(u))) ) * dΩ
end

function dc_mhd(x,dx,dy,α,dΩ)
  u, p, j, φ = x
  du , dp , dj , dφ  = dx
  v_u, v_p, v_j, v_φ = dy
  ∫( α*v_u⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u))) ) ) * dΩ
end

function ℓ_φ(dy,φ,n_Γ,dΓ)
  v_u, v_p, v_j, v_φ = dy
  ∫( -(v_j⋅n_Γ)*φ )*dΓ
end

function ℓ_f(dy,f,dΩ)
  v_u, v_p, v_j, v_φ = dy
  ∫( v_u⋅f )*dΩ
end

function ℓ_thin_wall(dy,cw,jw,n_Γ,dΓ)
  ∫( τ*(v_j⋅n_Γ)*jw ) * dΓ
end

function a_thin_wall(x,dy,cw,jw,n_Γ,dΓ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( τ*((v_j⋅n_Γ)*(j⋅n_Γ) + cw*(v_j⋅n_Γ)*(n_Γ⋅(∇(j)⋅n_Γ))) )*dΓ
end

function a_B(x,dy,γ,B,dΩ)
  u, p, j, φ = x
  v_u, v_p, v_j, v_φ = dy
  ∫( -(γ*(j×B)⋅v_u) - (u×B)⋅v_j )*dΩ
end


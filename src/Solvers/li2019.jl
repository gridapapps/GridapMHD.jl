
"""
  Implements a block-based solver for the MHD equations, 
  based on [(Li,2019)](https://doi.org/10.1137/19M1260372).

  There are three main components: 
    - A non-linear Newton solver, from GridapSolvers.jl 
    - A GMRES linear solver, from GridapSolvers.jl
    - A block-based upper-triangular preconditioner.
"""
function Li2019Solver(op::FEOperator,params)

  α1 = params[:fluid][:ζ] + params[:fluid][:β]
  inv_α1 = 1.0/α1

  # Preconditioner
  k  = max(params[:fespaces][:order_u],params[:fespaces][:order_j])
  model = params[:model]
  # Ωf = _interior(model,params[:fluid][:domain])
  # dΩ = Measure(Ωf,2*k)
  # a_j = li2019_Dj(dΩ,params)
  # a_u = li2019_Fk(dΩ,params)
  # a_p(p,v_p) = ∫(-inv_α1*p*v_p)*dΩ
  # a_φ(φ,v_φ) = ∫(-φ*v_φ)*dΩ
  a_j, a_u, a_p, a_φ = precond(params,k)

  U_u, U_p, U_j, U_φ = get_trial(op)
  V_u, V_p, V_j, V_φ = get_test(op)

  NB = 4
  diag_solvers = map(s -> get_block_solver(Val(s),params),params[:solver][:block_solvers])
  diag_blocks = [
    BiformBlock(a_j,U_j,V_j),
    TriformBlock(a_u,U_u,V_u,2),
    BiformBlock(a_p,U_p,V_p),
    BiformBlock(a_φ,U_φ,V_φ)
  ]
  blocks = map(CartesianIndices((NB,NB))) do I
    (I[1] == I[2]) ? diag_blocks[I[1]] : LinearSystemBlock()
  end
  coeffs = [1.0 2.0 0.0 2.0;  # | Dj  2Aju   0   2Ajφ |   | j |   | bj |
            0.0 1.0 1.0 0.0;  # | 0    Fk   Aup   0   |   | u |   | bu |
            0.0 0.0 1.0 0.0;  # | 0    0    Ip    0   | ⋅ | p | = | bp |
            0.0 0.0 0.0 1.0]  # | 0    0    0    Iφ   |   | φ |   | bφ |
  P = BlockTriangularSolver(blocks,diag_solvers,coeffs,:upper)

  # Linear Solver
  verbose = i_am_main(get_parts(model))
  nl_rtol = params[:solver][:rtol]
  l_rtol  = nl_rtol/10.0
  
  m = params[:solver][:niter_ls]
  l_solver = FGMRESSolver(m,P;rtol=l_rtol,atol=1e-14,verbose=verbose)

  # Nonlinear Solver
  niter = params[:solver][:niter]
  nl_solver = GridapSolvers.NewtonSolver(l_solver,maxiter=niter,atol=1e-14,rtol=nl_rtol,verbose=verbose)
  return nl_solver
end

function precond(params,k)
  dΩf, α, β, γ, σf, f, B, ζ, g = retrieve_fluid_params(params)
  solid = params[:solid]
  dΩs, σs = retrieve_solid_params(params)

  function a_u(u,du,dv)
    r = a_mhd_u_u(du,dv,β,dΩf) + n_dc_mhd_u_u(u,du,dv,α,dΩf) + ∫(γ⋅(du×B)⋅(dv×B)) * dΩf
    if abs(ζ) > eps(typeof(ζ))
      r = r + ∫( ζ*(∇⋅du)*(∇⋅dv) ) * dΩf
    end
    return r
  end

  function a_j(j,v_j) 
    r = ∫( (j⋅v_j) + (∇⋅j)⋅(∇⋅v_j))*dΩf 
    if solid !== nothing
      r = r + ∫( (j⋅v_j) + (∇⋅j)⋅(∇⋅v_j))*dΩs
    end
    return r
  end

  α1 = ζ+β
  inv_α1 = 1.0/α1
  a_p(p,v_p) = ∫(-inv_α1*p*v_p)*dΩf

  function a_φ(φ,v_φ)
    r =  ∫(-φ*v_φ)*dΩf
    if solid !== nothing
      r = r +∫(-φ*v_φ)*dΩs
    end
    return r
  end

  return a_j, a_u, a_p, a_φ
end

function li2019_Dj(dΩ,params)
  k = params[:fespaces][:k]
  _, params_thin_wall, _, _ = retrieve_bcs_params(params,k)
  function a_j(j,v_j) 
    r = a_mhd_j_j(j,v_j,dΩ) + ∫((∇⋅j)⋅(∇⋅v_j))*dΩ 
    for p in params_thin_wall
      τ,cw,jw,n_Γ,dΓ = p
      r += a_thin_wall_j_j(j,v_j,τ,cw,jw,n_Γ,dΓ)
    end
    return r
  end
  return a_j
end

function li2019_Fk(dΩ,params)
  fluid = params[:fluid]
  α, β, γ, B  = fluid[:α], fluid[:β], fluid[:γ], fluid[:B]
  a_fk(u,du,dv) = a_mhd_u_u(du,dv,β,dΩ) + n_dc_mhd_u_u(u,du,dv,α,dΩ) + ∫(γ⋅(du×B)⋅(dv×B)) * dΩ
  return a_fk
end

function li2019_laplacian(dΩ,params)
  p_conformity = params[:fespaces][:p_conformity]

  model = params[:model]
  Γ  = Boundary(model)
  Λ  = Skeleton(model)

  k  = params[:fespaces][:k]
  dΓ = Measure(Γ,2*k)
  dΛ = Measure(Λ,2*k)

  n_Γ = get_normal_vector(Γ)
  n_Λ = get_normal_vector(Λ)

  h_e_Λ = get_cell_size(Λ)
  h_e_Γ = get_cell_size(Γ)

  β = 100.0
  function a_Δp(u,v)
    r = ∫(∇(u)⋅∇(v))*dΩ
    if p_conformity == :L2
      r += ∫(-jump(u⋅n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ) + β/h_e_Λ*jump(u⋅n_Λ)⋅jump(v⋅n_Λ))*dΛ
      r += ∫(-(∇(u)⋅n_Γ)⋅v - u⋅(∇(v)⋅n_Γ) + β/h_e_Γ*(u⋅n_Γ)⋅(v⋅n_Γ))*dΓ
    end
    return r
  end
  return a_Δp
end

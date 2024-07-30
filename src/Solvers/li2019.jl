
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
  model = params[:model]
  k  = max(params[:fespaces][:order_u],params[:fespaces][:order_j])
  Ωf = _interior(model,params[:fluid][:domain])
  dΩ = Measure(Ωf,2*k)
  a_Dj = li2019_Dj(dΩ,params)
  a_Fk = li2019_Fk(dΩ,params)
  a_Ip(p,v_p) = ∫(-inv_α1*p*v_p)*dΩ
  a_Iφ(φ,v_φ) = ∫(-φ*v_φ)*dΩ

  U_u, U_p, U_j, U_φ = get_trial(op)
  V_u, V_p, V_j, V_φ = get_test(op)

  diag_solvers = map(s -> get_block_solver(Val(s),params),params[:solver][:block_solvers])

  j_block = BiformBlock(a_Dj,U_j,V_j)
  u_block = TriformBlock(a_Fk,U_u,U_u,V_u,[2])
  p_block = BiformBlock(a_Ip,U_p,V_p)
  φ_block = BiformBlock(a_Iφ,U_φ,V_φ)
  blocks = [    j_block        LinearSystemBlock() LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock()     u_block         LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock()      p_block         LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock() LinearSystemBlock()      φ_block       ]
  coeffs = [1.0 2.0 0.0 2.0;  # | Dj  2Aju   0   2Aφu |   | j |   | bj |
            0.0 1.0 1.0 0.0;  # | 0    Fk   Aup   0   |   | u |   | bu |
            0.0 0.0 1.0 0.0;  # | 0    0    Ip    0   | ⋅ | p | = | bp |
            0.0 0.0 0.0 1.0]  # | 0    0    0    Iφ   |   | φ |   | bφ |
  P = BlockTriangularSolver(blocks,diag_solvers,coeffs,:upper)

  # Linear Solver
  verbose = i_am_main(get_parts(model))
  nl_rtol = params[:solver][:rtol]
  l_rtol  = nl_rtol/10.0
  
  m = params[:solver][:niter]
  l_solver = FGMRESSolver(m,P;rtol=l_rtol,atol=1e-14,verbose=verbose)

  # Nonlinear Solver
  nl_solver = GridapSolvers.NewtonSolver(l_solver,maxiter=10,atol=1e-14,rtol=nl_rtol,verbose=verbose)
  return nl_solver
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
  a_fk(u,du,dv) = a_mhd_u_u(du,dv,β,dΩ) + dc_mhd_u_u(u,du,dv,α,dΩ) + ∫(γ⋅(du×B)⋅(dv×B)) * dΩ
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

  h_e_Λ = get_edge_measures(Λ,dΛ)
  h_e_Γ = get_edge_measures(Γ,dΓ)

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

get_edge_measures(Ω::Triangulation,dΩ) = sqrt∘CellField(get_array(∫(1)dΩ),Ω)
function get_edge_measures(Ω::GridapDistributed.DistributedTriangulation,dΩ)
  return sqrt∘CellField(map(get_array,local_views(∫(1)*dΩ)),Ω)
end

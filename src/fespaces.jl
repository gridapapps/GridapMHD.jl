
# MultiFieldStyle

_multi_field_style(params) = _multi_field_style(Val(params[:solver][:solver]))
_multi_field_style(::Val{:julia}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:petsc}) = ConsecutiveMultiFieldStyle()
_multi_field_style(::Val{:li2019}) = BlockMultiFieldStyle(4,(1,1,1,1),(3,1,2,4)) # (j,u,p,φ)
_multi_field_style(::Val{:badia2024}) = BlockMultiFieldStyle(3,(2,1,1),(1,3,2,4)) # ([u,j],p,φ)
_multi_field_style(::Val{:h1h1block}) = BlockMultiFieldStyle(3,(1,1,1),(1,2,3)) # (u,p,φ)

# FE Spaces

function setup_fe_spaces(params)
  formulation = params[:fespaces][:formulation]
  
  if formulation ∈ (:H1HDiv,:HDivHDiv)
    U_u, V_u = fe_space_u(params)
    U_p, V_p = fe_space_p(params)
    U_j, V_j = fe_space_j(params)
    U_φ, V_φ = fe_space_φ(params)

    mfs = _multi_field_style(params)
    V = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=mfs)
    if !has_transient(params)
      U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)
    else
      U = TransientMultiFieldFESpace([U_u,U_p,U_j,U_φ];style=mfs)
    end
  elseif formulation ∈ (:H1H1,)
    U_u, V_u = fe_space_u(params)
    U_p, V_p = fe_space_p(params)
    U_φ, V_φ = fe_space_φ(params)

    mfs = _multi_field_style(params)
    V = MultiFieldFESpace([V_u,V_p,V_φ];style=mfs)
    if !has_transient(params)
      U = MultiFieldFESpace([U_u,U_p,U_φ];style=mfs)
    else
      U = TransientMultiFieldFESpace([U_u,U_p,U_φ];style=mfs)
    end
  else
    @error "Unknown formulation: $formulation"
  end

  return U, V
end

function fe_space_u(params)
  uses_mg = space_uses_multigrid(params[:solver])[1]

  Ωf = uses_mg ? params[:multigrid][:Ωf] : params[:Ωf]

  reffe_u = params[:fespaces][:reffe_u]
  u_bc = params[:bcs][:u][:values]
  V_u = TestFESpace(Ωf,reffe_u;dirichlet_tags=params[:bcs][:u][:tags])
  U_u = _trial_fe_space(V_u,u_bc,params)

  if uses_mg
    params[:multigrid][:trials][:u] = U_u
    params[:multigrid][:tests][:u] = V_u
    U_u, V_u = get_fe_space(U_u,1), get_fe_space(V_u,1)
  end

  return U_u, V_u
end

function fe_space_p(params)
  @notimplementedif space_uses_multigrid(params[:solver])[2]

  Ωf = params[:Ωf]
  conformity = params[:fespaces][:p_conformity]
  constraint = params[:fespaces][:p_constraint]

  reffe_p = params[:fespaces][:reffe_p]
  V_p = TestFESpace(Ωf,reffe_p;conformity,constraint)
  U_p = _trial_fe_space(V_p,nothing,params)

  return U_p, V_p
end

function fe_space_j(params)
  uses_mg = space_uses_multigrid(params[:solver])[3]
  model = uses_mg ? params[:multigrid][:mh] : params[:model]

  reffe_j = params[:fespaces][:reffe_j]
  j_bc = params[:bcs][:j][:values]
  V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
  U_j = _trial_fe_space(V_j,j_bc,params)

  if uses_mg
    params[:multigrid][:trials][:j] = U_j
    params[:multigrid][:tests][:j] = V_j
    U_j, V_j = get_fe_space(U_j,1), get_fe_space(V_j,1)
  end

  return U_j, V_j
end

function fe_space_φ(params)
  @notimplementedif space_uses_multigrid(params[:solver])[4]
  model = params[:model]

  reffe_φ = params[:fespaces][:reffe_φ]
  conformity = params[:fespaces][:φ_conformity]
  constraint = params[:fespaces][:φ_constraint]
  if conformity != :H1
    dirichlet_tags = params[:bcs][:φ][:tags]
    φ_bc = params[:bcs][:φ][:values]
    V_φ = TestFESpace(model,reffe_φ;conformity,constraint)
    U_φ = _trial_fe_space(V_φ,nothing,params)
  else
    dirichlet_tags = params[:bcs][:φ][:tags]
    φ_bc = params[:bcs][:φ][:values]
    V_φ = TestFESpace(model,reffe_φ;conformity,constraint,dirichlet_tags)
    U_φ = _trial_fe_space(V_φ,φ_bc,params)
  end

  return U_φ, V_φ
end

function _trial_fe_space(V,v_bc,params)
  if isnothing(v_bc) || (isa(v_bc,Number) && iszero(v_bc))
    return V
  elseif has_transient(params)
    return TransientTrialFESpace(V,v_bc)
  else
    return TrialFESpace(V,v_bc)
  end
end

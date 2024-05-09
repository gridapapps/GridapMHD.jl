
# Generic methods for parameter checking

function _check_mandatory(params,mandatory,p)
  if !isa(params,Dict{Symbol})
    error("The params$p has to be a Dict{Symbol}")
  end
  for key in keys(mandatory)
    if mandatory[key] && !haskey(params,key)
      error("Key :$key is a mandatory key in params$p, but it is not provided.")
    end
  end
end

function _add_optional(_subparams,mandatory,optional,params,p)
  # New dict
  subparams = Dict{Symbol,Any}()
  merge!(subparams,_subparams)
  # Check that we have computed defaults for all optionals
  for key in keys(mandatory)
    if !mandatory[key] && !haskey(optional,key)
      error("Internal error")
    end
  end
  # Set default args
  for key in keys(optional)
    if !haskey(subparams,key)
      subparams[key] = optional[key]
    end
  end
  subparams
end

function _add_optional_weak(_subparams,optional,params,p)
  # New dict
  subparams = Dict{Symbol,Any}()
  merge!(subparams,_subparams)
  # Set default args
  for key in keys(optional)
    if !haskey(subparams,key)
      subparams[key] = optional[key]
    end
  end
  subparams
end

function _check_unused(subparams,mandatory,params,p)
  if params[:check_valid]
    for key in keys(subparams)
      if !haskey(mandatory,key)
        error("Key :$key is not a valid key in params$p. Set params[:check_valid] = false to ignore invalid keys.")
      end
    end
  end
end

function _check_mandatory_and_add_optional(_subparams,mandatory,optional,params,p)
  _check_mandatory(_subparams,mandatory,p)
  subparams = _add_optional(_subparams,mandatory,optional,params,p)
  _check_unused(subparams,mandatory,params,p)
  subparams
end

function _check_mandatory_and_add_optional_weak(t,mandatory,optional,params,p)
  function _params_bcs(_t)
    error("The value params$p has to be a Dict{Symbol} or a vector of Dict{Symbol}.")
  end
  function _params_bcs(_t::AbstractVector)
    map(i->_params_bcs(_t[i],"[$i]"),1:length(_t))
  end
  function _params_bcs(_t::Dict{Symbol},i="")
    _check_mandatory_and_add_optional(_t,mandatory,optional,params,"$p$i")
  end
  r = _params_bcs(t)
  isa(r,Dict) ? [r] : r
end

"""
    add_default_params(params::Dict{Symbol}) -> full_params

Return a new `Dict` with the contents of `params` plus
default values for the optional keys not provided in `params`.
It also checks the validity of the main parameter dictionary `params`.
"""
function add_default_params(_params)
  mandatory = Dict(
    :ptimer=>false,
    :debug=>false,
    :solve=>true,
    :res_assemble=>false,
    :jac_assemble=>false,
    :model=>true,
    :fluid=>true,
    :solid=>false,
    :bcs=>true,
    :fespaces=>false,
    :solver=>true,
    :multigrid=>false,
    :check_valid=>false,
  )
  _check_mandatory(_params,mandatory,"")
  optional = Dict(
    :ptimer=>default_ptimer(_params[:model]),
    :debug=>false,
    :res_assemble=>false,
    :jac_assemble=>false,
    :solid=>nothing,
    :fespaces=>nothing,
    :multigrid=>nothing,
    :check_valid=>true,
  )
  params = _add_optional(_params,mandatory,optional,_params,"")
  _check_unused(params,mandatory,params,"")
  # Process sub-params
  params[:fluid] = params_fluid(params)
  if params[:solid] !== nothing
    params[:solid] = params_solid(params)
  end
  params[:bcs] = params_bcs(params)
  params[:fespaces] = params_fespaces(params)
  params[:solver] = params_solver(params)
  #params[:multigrid] = params_multigrid(params)
  params
end

default_ptimer(model) = PTimer(DebugArray(LinearIndices((1,))))
default_ptimer(model::GridapDistributed.DistributedDiscreteModel) = PTimer(get_parts(model))

"""
Valid keys for `params[:solver]` are the following:

# Mandatory keys
- `:solver`: Name of the selected solver.
            Valid values are `[:julia, :petsc, :li2019, :badia2024]`.

# Optional keys
-  `:rtol` : Relative tolerance.
-  `:matrix_type`: Matrix type for the linear system.
-  `:vector_type`: Vector type for the linear system.
-  `:solver_postpro`: Function to postprocess the solver cache, with signature f(cache,info)

# Solver-dependent optional keys
-  `:niter`: Number of iterations for the linear solver (only for iterative solvers).
-  `:petsc_options`: PETSc options for the linear solver (only if PETSc is used).
-  `:block_solvers`: Array of solvers for the diagonal blocks (only for block-based solvers).
"""
function params_solver(params::Dict{Symbol,Any})
  if isa(params[:solver],Symbol)
    solver = default_solver_params(Val(params[:solver]))
    return solver
  end
  @assert haskey(params[:solver],:solver)

  optional = default_solver_params(Val(params[:solver][:solver]))
  solver   = _add_optional_weak(params[:solver],optional,params,"[:solver]")
  return solver
end

function default_solver_params(::Val{:julia})
  return Dict(
    :solver => :julia,
    :matrix_type    => SparseMatrixCSC{Float64,Int64},
    :vector_type    => Vector{Float64},
    :solver_postpro => ((cache,info) -> nothing),
    :rtol           => 1e-5,
  )
end

function default_solver_params(::Val{:petsc})
  return Dict(
    :solver => :petsc,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :solver_postpro => ((cache,info) -> snes_postpro(cache,info)),
    :petsc_options  => "-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_7 0",
    :niter          => 100,
    :rtol           => 1e-5,
  )
end

function default_solver_params(::Val{:li2019})
  return Dict(
    :solver => :li2019,
    :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
    :vector_type    => Vector{PetscScalar},
    :petsc_options  => "-ksp_error_if_not_converged true -ksp_converged_reason",
    :solver_postpro => ((cache,info) -> gridap_postpro(cache,info)),
    :block_solvers  => [:petsc_mumps,:petsc_gmres_schwarz,:petsc_cg_jacobi,:petsc_cg_jacobi],
    :niter          => 80,
    :rtol           => 1e-5,
  )
end

uses_petsc(solver::Dict) = uses_petsc(solver[:solver])
uses_petsc(solver::Symbol) = uses_petsc(Val(solver))
uses_petsc(::Val{:julia}) = false
uses_petsc(::Val{:petsc}) = true
uses_petsc(::Val{:li2019}) = true

uses_multigrid(solver::Dict) = any(space_uses_multigrid(solver))
space_uses_multigrid(solver::Dict) = space_uses_multigrid(Val(solver[:solver]),solver)
space_uses_multigrid(::Val{:julia},solver) = fill(false,4)
space_uses_multigrid(::Val{:petsc},solver) = fill(false,4)
space_uses_multigrid(::Val{:li2019},solver) = map(s -> s==:gmg, solver[:block_solvers])

function snes_postpro(cache,info)
  snes = cache.snes[]
  i_petsc = Ref{PetscInt}()
  @check_error_code GridapPETSc.PETSC.SNESGetIterationNumber(snes,i_petsc)
  info[:nls_iters] = Int(i_petsc[])
  @check_error_code GridapPETSc.PETSC.SNESGetLinearSolveIterations(snes,i_petsc)
  info[:ls_iters] = Int(i_petsc[])
  nothing
end

function gridap_postpro(cache,info)
  ls = cache.ns.solver
  log = ls.log

  info[:ls_iters] = log.num_iters
  info[:ls_residuals] = log.residuals[1:log.num_iters+1]
end

"""
Valid keys for `params[:fespaces]` are the following:

# Optional keys
-  `:k`: Polynomial degree for the fluid velocity.
-  `:p_space`: FESpace conformity for pressure. Possible values are [:P,:Q]
"""
function params_fespaces(params::Dict{Symbol,Any})
  if !haskey(params,:fespaces) || isa(params[:fespaces],Nothing)
    params[:fespaces] = Dict{Symbol,Any}()
  end
  mandatory = Dict(
   :k => false,
   :p_space => false,
  )
  optional = Dict(
   :k => 2,
   :p_space => :P,
  )
  fespaces = _add_optional(params[:fespaces],mandatory,optional,params,"[:fespaces]")
  fespaces[:p_conformity] = p_conformity(params[:model],fespaces)
  return fespaces
end

function p_conformity(poly::Polytope,feparams)
  p_space = feparams[:p_space]
  if is_n_cube(poly) && (p_space == :P)
    conf = :L2
  else
    conf = :H1
  end
  return conf
end
function p_conformity(Ω::Triangulation,feparams)
  reffes = get_reffes(Ω)
  @assert length(reffes) == 1
  return p_conformity(get_polytope(first(reffes)),feparams)
end
function p_conformity(Ω::GridapDistributed.DistributedTriangulation,feparams)
  p = map(local_views(Ω)) do Ωi
    p_conformity(Ωi,feparams)
  end
  return getany(p) # We assume same polytope in all parts
end
p_conformity(model::DiscreteModel,feparams) = p_conformity(Interior(model),feparams)
p_conformity(model::GridapDistributed.DistributedDiscreteModel,feparams) = p_conformity(Interior(model),feparams)

"""
Valid keys for `params[:multigrid]` are the following:

# Optional keys
-  `:k`: Polynomial degree for the fluid velocity.
-  `:p_space`: FESpace conformity for pressure. Possible values are [:P,:Q]
"""
function params_multigrid(params::Dict{Symbol,Any})
  solver = params[:solver]
  if !uses_multigrid(solver)
    if !isa(params[:multigrid],Nothing)
      @warn "Multigrid is not used with solver $(solver[:solver]). Ignoring params[:multigrid]."
    end
    return nothing
  end
  if isa(params[:multigrid],Nothing) || !haskey(params[:multigrid],:mh)
    @error "Multigrid is used with solver $(solver[:solver]), but params[:multigrid] is missing!"
  end
  multigrid = params[:multigrid]
  return multigrid
end

"""
Valid keys for `params[:fluid]` are the following.

# Mandatory keys
- `domain`: Domain where to solve the MHD problem.
  The domain is represented with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
  It can also be a `Gridap.DiscreteModel` or `GridapDistributed.DistributedDiscreteModel`
  if `params[:solid]` is `nothing`. In this case,
  one needs to guarantee `params[:model]===params[:fluid][:domain]`.
-  `:α`: Value of the parameter `α`.
-  `:β`: Value of the parameter `β`.
-  `:γ`: Value of the parameter `γ`.
-  `:B`: Value of the parameter `B`.

# Optional keys
-  `:f=>VectorValue(0,0,0)`: Value of the parameter `f`.
-  `:σ=>1`: Value of the parameter `σ`.
-  `:ζ=>0`: Value of the augmented lagrangian weight.
"""
function params_fluid(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :α=>true,
   :β=>true,
   :γ=>true,
   :B=>true,
   :f=>false,
   :σ=>false,
   :ζ=>false,
  )
  optional = Dict(:σ=>1.0,:f=>VectorValue(0,0,0),:ζ=>0.0)
  fluid = _check_mandatory_and_add_optional(params[:fluid],mandatory,optional,params,"[:fluid]")
  fluid
end

"""
Valid keys for `params[:solid]` are the following

# Mandatory keys
- `domain`: Domain occupied by the solid.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
# Optional keys
-  `:σ=>1`: Value of the parameter `σ`.
"""
function params_solid(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :σ=>false,
  )
  optional = Dict(:σ=>1)
  solid = _check_mandatory_and_add_optional(params[:solid],mandatory,optional,params,"[:solid]")
  solid
end

"""
Valid keys for `params[:bcs]` are the following

# Mandatory keys
- `:u`: A `Dict` defining strong Dirichlet conditions on the fluid velocity.
   See  [`params_bcs_u`](@ref) for further details.
- `:j`: A `Dict` defining strong Dirichlet conditions on the charge current.
   See  [`params_bcs_j`](@ref) for further details.

# Optional keys
- `:φ=>[]`: A `Dict` or a vector of `Dict`s.
   Each `Dict` defines a weak boundary condition for the potential.
   See  [`params_bcs_φ`](@ref) for further details.
- `:t=>[]`: A `Dict` or a vector of `Dict`s.
  Each `Dict` defines a boundary traction on the fluid.
   See  [`params_bcs_t`](@ref) for further details.
- `:thin_wall=>[]`: A `Dict` or a vector of `Dict`s.
  Each Dict defines a thin wall law.
   See  [`params_bcs_thin_wall`](@ref) for further details.
"""
function params_bcs(params)
  mandatory = Dict(
   :u=>true,
   :j=>true,
   :φ=>false,
   :t=>false,
   :thin_wall=>false,
   :f => false,
   :B => false,
  )
  optional = Dict(
   :φ=>[],
   :t=>[],
   :thin_wall=>[],
   :f =>[],
   :B =>[],
  )
  bcs = _check_mandatory_and_add_optional(params[:bcs],mandatory,optional,params,"[:bcs]")
  # Sub params
  bcs[:u] = params_bcs_u(params)
  bcs[:j] = params_bcs_j(params)
  if bcs[:φ] !== optional[:φ]
    bcs[:φ] = params_bcs_φ(params)
  end
  if bcs[:t] !== optional[:t]
    bcs[:t] = params_bcs_t(params)
  end
  if bcs[:thin_wall] !== optional[:thin_wall]
    bcs[:thin_wall] = params_bcs_thin_wall(params)
  end
  bcs
end

"""
Valid keys for `params[:bcs][:u]` are the following

# Mandatory keys
- `:tags`: Dirichlet tags where to impose strong boundary conditions for the fluid velocity.

# Optional keys
- `:values => zero_values(prams[:bcs][:u][:tags])`: The fluid velocity value or function
   to be imposed at each of the given tags.
"""
function params_bcs_u(params::Dict{Symbol,Any})
  mandatory = Dict(
   :tags=>true,
   :values=>false,
  )
  _check_mandatory(params[:bcs][:u],mandatory,"[:bcs][:u]")
  optional = Dict(
    :values=>zero_values(params[:bcs][:u][:tags]),
  )
  u = _add_optional(params[:bcs][:u],mandatory,optional,params,"[:bcs][:u]")
  _check_unused(u,mandatory,params,"[:bcs][:u]")
  u
end

zero_values(tags) = VectorValue(0,0,0)
zero_values(tags::Vector) = map(zero_values,tags)

"""
Valid keys for `params[:bcs][:j]` are the following

# Mandatory keys
- `:tags`: Dirichlet tags where to impose strong boundary conditions for the charge current in normal direction.

# Optional keys
- `:values => zero_values(params[:bcs][:j][:tags])`: The charge current value or function
   to be imposed at each of the given tags.
"""
function params_bcs_j(params::Dict{Symbol,Any})
  mandatory = Dict(
   :tags=>true,
   :values=>false,
  )
  _check_mandatory(params[:bcs][:j],mandatory,"[:bcs][:j]")
  optional = Dict(
    :values=>zero_values(params[:bcs][:j][:tags]),
  )
  j = _add_optional(params[:bcs][:j],mandatory,optional,params,"[:bcs][:j]")
  _check_unused(j,mandatory,params,"[:bcs][:j]")
  j
end

"""
Valid keys for the dictionaries in `params[:bcs][:φ]` are the following.

# Mandatory keys
- `:domain`: Domain where to impose the potential weakly.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:value`: Value of the electric potential to be imposed weakly.
"""
function params_bcs_φ(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :value=>true,
  )
  optional = Dict()
  _check_mandatory_and_add_optional_weak(params[:bcs][:φ],mandatory,optional,params,"[:bcs][:φ]")
end

"""
Valid keys for the dictionaries in `params[:bcs][:t]` are the following.

# Mandatory keys
- `:domain`: Domain where to impose the fluid boundary traction weakly.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:value`: Value of the fluid traction to be imposed weakly.
"""
function params_bcs_t(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :value=>true,
  )
  optional = Dict()
  _check_mandatory_and_add_optional_weak(params[:bcs][:t],mandatory,optional,params,"[:bcs][:t]")
end

"""
Valid keys for the dictionaries in `params[:bcs][:thin_wall]` are the following.

The thin wall law is

    j⋅n + cw*n⋅∇(j)⋅n = jw

where `j` is unknown, `n`  is the boundary outward normal, and `cw,jw` are parameters.
The thin wall law is imposed weakly via a penalty parameter `τ`.

# Mandatory keys
- `:domain`: Domain where to impose the thin wall law.
  The domain is represented either with a `Triangulation` or with a `Integer`/`String`
  tag in the underlying discrete model.
- `:cw`: Value of the parameter `cw`.
- `:τ`: Value of the parameter `τ`.

# Optional keys
- `:jw=>0`: Value of the parameter `jw`.
"""
function params_bcs_thin_wall(params::Dict{Symbol,Any})
  mandatory = Dict(
   :domain=>true,
   :cw=>true,
   :τ=>true,
   :jw=>false,
  )
  optional = Dict(:jw=>0)
  _check_mandatory_and_add_optional_weak(params[:bcs][:thin_wall],mandatory,optional,params,"[:bcs][:thin_wall]")
end

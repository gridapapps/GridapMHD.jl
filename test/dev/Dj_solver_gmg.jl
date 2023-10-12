
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
using GridapSolvers.MultilevelTools
import GridapSolvers.PatchBasedSmoothers as PBS
using Gridap.ReferenceFEs, Gridap.Geometry

using GridapMHD

function setup(
  model;
  N  = 1.0,
  Ha = 1.0,
  cw = 0.0,
  τ  = 100
  )

  α = (1.0/N)
  β = (1.0/Ha^2)
  γ = 1.0

  u_inlet((x,y,z)) = VectorValue(36.0*(y-1/4)*(y+1/4)*(z-1)*(z+1),0,0)

  _params = Dict(
    :debug=>false,
    :solve=>true,
    :res_assemble=>false,
    :jac_assemble=>false,
    :model => model,
    :fluid=>Dict(
      :domain=>model,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :f=>VectorValue(0.0,0.0,0.0),
      :B=>VectorValue(0.0,1.0,0.0),
    ),
    :solver=>:julia,
   )

  if cw == 0.0
   _params[:bcs] = Dict( 
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet, VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
		    :tags => ["wall", "inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)],
      )
    )

  else 
   _params[:bcs] = Dict(
      :u => Dict(
        :tags => ["inlet", "wall"],
        :values => [u_inlet, VectorValue(0.0, 0.0, 0.0)]
      ),
      :j => Dict(
        :tags => ["inlet", "outlet"], 
        :values=>[VectorValue(0.0,0.0,0.0), VectorValue(0.0,0.0,0.0)]
      ),
      :thin_wall => [Dict(
        :τ=>τ,
        :cw=>cw,
        :domain => Boundary(model, tags="wall")
      )]
    )
  end

  params = GridapMHD.add_default_params(_params)
  return params
end

function _Dj(dΩ,params)
  fluid = params[:fluid]
  γ  = fluid[:γ]
  k = params[:fespaces][:k]

  params_thin_wall = []
  bcs = params[:bcs]
  for i in 1:length(bcs[:thin_wall])
    τ   = bcs[:thin_wall][i][:τ]
    cw  = bcs[:thin_wall][i][:cw]
    jw  = bcs[:thin_wall][i][:jw]
    # Γ   = BoundaryTriangulation(PD)
    #dΓ  = Measure(Γ,2*k)
    #n_Γ = get_normal_vector(Γ)
    #push!(params_thin_wall,(τ,cw,jw,n_Γ,dΓ))
  end

  function a_j(j,v_j) 
    r = ∫(j⋅v_j + (∇⋅j)⋅(∇⋅v_j))*dΩ 
    #for p in params_thin_wall
    #  τ,cw,jw,n_Γ,dΓ = p
    #  r += ∫(τ*(v_j⋅n_Γ)⋅(j⋅n_Γ) + cw*(v_j⋅n_Γ)⋅(n_Γ⋅(∇(j)⋅n_Γ)))*dΓ
    #end
    return r
  end
  return a_j
end

function test_solver(s,D_j)
  ns = numerical_setup(symbolic_setup(s,D_j),D_j)

  b = GridapSolvers.allocate_col_vector(D_j)
  x = GridapSolvers.allocate_col_vector(D_j)

  fill!(b,1.0)
  solve!(x,ns,b)
  err = norm(b - D_j*x)

  return err
end

function test_smoother(s,D_j)
  ns = numerical_setup(symbolic_setup(s,D_j),D_j)
  b = GridapSolvers.allocate_col_vector(D_j)
  x = GridapSolvers.allocate_col_vector(D_j)
  r = GridapSolvers.allocate_row_vector(D_j)
  fill!(b,1.0)
  fill!(x,1.0)
  mul!(r,D_j,x)
  r .= b .- r
  solve!(x,ns,r)
  err = norm(b - D_j*x)
  return err
end

function PBS.PatchDecomposition(mh::ModelHierarchy)
  pds = Vector{PBS.DistributedPatchDecomposition}(undef,num_levels(mh))
  for lev in 1:num_levels(mh)
    model = get_model(mh,lev)
    pds[lev] = PBS.PatchDecomposition(model)
  end
  return pds
end

function get_hierarchy_matrices(mh,tests,trials)
  mats = Vector{AbstractMatrix}(undef,num_levels(mh))
  A = nothing
  b = nothing
  for lev in 1:num_levels(mh)
    model = get_model(mh,lev)
    U_j = get_fe_space(trials,lev)
    V_j = get_fe_space(tests,lev)
    Ω   = Triangulation(model)
    dΩ  = Measure(Ω,2*k)
    ai  = GridapMHD.Li2019._Dj(dΩ,params)
    if lev == 1
      f = VectorValue(1.0,1.0,1.0)
      li(v) = ∫(v⋅f)*dΩ
      op    = AffineFEOperator(ai,li,U_j,V_j)
      A, b  = get_matrix(op), get_vector(op)
      mats[lev] = A
    else
      mats[lev] = assemble_matrix(ai,U_j,V_j)
    end
  end
  return mats, A, b
end

function get_patch_smoothers(tests,patch_spaces,patch_decompositions,qdegree,params)
  mh = tests.mh
  nlevs = num_levels(mh)
  smoothers = Vector{RichardsonSmoother}(undef,nlevs-1)
  for lev in 1:nlevs-1
    parts = get_level_parts(mh,lev)
    if i_am_in(parts)
      PD = patch_decompositions[lev]
      Ph = get_fe_space(patch_spaces,lev)
      Vh = get_fe_space(tests,lev)
      Ω  = Triangulation(PD)
      dΩ = Measure(Ω,qdegree)
      a_j = _Dj(dΩ,params)
      local_solver = LUSolver() # IS_ConjugateGradientSolver(;reltol=1.e-6)
      patch_smoother = PatchBasedLinearSolver(a_j,Ph,Vh,local_solver)
      smoothers[lev] = RichardsonSmoother(patch_smoother,10,1.0)
    end
  end
  return smoothers
end

############################################################################################

np = 1
ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

# Geometry
#model  = GridapMHD.Meshers.expansion_generate_base_mesh()
mh = GridapMHD.Meshers.expansion_generate_mesh_hierarchy(ranks,1,[1,1,1]);
model = get_model(mh,1)
params = setup(model);

# FESpaces
k = 1#params[:fespaces][:k]
qdegree = 2*k+2
j_bc = params[:bcs][:j][:values]
reffe_j = ReferenceFE(raviart_thomas,Float64,k-1)
tests  = TestFESpace(mh,reffe_j;dirichlet_tags=params[:bcs][:j][:tags]);
trials = TrialFESpace(tests,j_bc);

# Patch solver
patch_decompositions = PBS.PatchDecomposition(mh)
patch_spaces = PBS.PatchFESpace(mh,reffe_j,DivConformity(),patch_decompositions,tests);
smoothers = get_patch_smoothers(tests,patch_spaces,patch_decompositions,qdegree,params)

restrictions, prolongations = setup_transfer_operators(trials,qdegree;mode=:residual);
smatrices, A, b = get_hierarchy_matrices(mh,tests,trials);

gmg = GMGLinearSolver(mh,
                      smatrices,
                      prolongations,
                      restrictions,
                      pre_smoothers=smoothers,
                      post_smoothers=smoothers,
                      maxiter=4,
                      rtol=1.0e-8,
                      verbose=true,
                      mode=:preconditioner)

solver = FGMRESSolver(100,gmg;rtol=1e-6,verbose=true)

ns = numerical_setup(symbolic_setup(solver,A),A)
x = GridapSolvers.allocate_col_vector(A)
solve!(x,ns,b)


Pl = LinearSolvers.IdentitySolver()
solver2 = GMRESSolver(1000;Pl=Pl,rtol=1e-6,verbose=true)
ns2 = numerical_setup(symbolic_setup(solver2,A),A)
x2 = GridapSolvers.allocate_col_vector(A)
solve!(x2,ns2,b)

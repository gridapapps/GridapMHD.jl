
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
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

function _Dj(PD,dΩ,params)
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
    for p in params_thin_wall
      τ,cw,jw,n_Γ,dΓ = p
      r += ∫(τ*(v_j⋅n_Γ)⋅(j⋅n_Γ) + cw*(v_j⋅n_Γ)⋅(n_Γ⋅(∇(j)⋅n_Γ)))*dΓ
    end
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

############################################################################################

np = 1
ranks = with_mpi() do distribute
  distribute(LinearIndices((np,)))
end

# Geometry
#model  = GridapMHD.Meshers.expansion_generate_base_mesh()
mh = GridapMHD.Meshers.expansion_generate_mesh_hierarchy(ranks,2,[1,1])
params = setup(model);

# FESpaces
k = params[:fespaces][:k]
j_bc = params[:bcs][:j][:values]
reffe_j = ReferenceFE(raviart_thomas,Float64,k-1)
V_j = TestFESpace(model,reffe_j;dirichlet_tags=params[:bcs][:j][:tags])
U_j = TrialFESpace(V_j,j_bc)

Ω = Triangulation(model)
dΩ = Measure(Ω,2*k)
a_j = GridapMHD.Li2019._Dj(dΩ,params)
D_j = assemble_matrix(a_j,U_j,V_j)

# Patch solver
PD  = PBS.PatchDecomposition(model)
P_j = PBS.PatchFESpace(model,reffe_j,DivConformity(),PD,V_j)

Ωp   = Triangulation(PD)
dΩp  = Measure(Ωp,2*k)
ap_j = _Dj(PD,dΩp,params)

local_solver = LUSolver()
patch_solver = PatchBasedLinearSolver(ap_j,P_j,U_j,local_solver)

#smoother = RichardsonSmoother(patch_solver,10,1.0)
#test_smoother(smoother,D_j)

solver = FGMRESSolver(100,patch_solver;rtol=1e-6,verbose=true)
test_solver(solver,D_j)

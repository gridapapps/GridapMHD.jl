using FileIO
using BSON
using Gridap, GridapDistributed, GridapSolvers
using Gridap.MultiField
using Gridap.Algebra
using PartitionedArrays

using GridapMHD: weak_form, add_default_params

t_parts = get_part_ids(SequentialBackend(),1)
t = PTimer(t_parts,verbose=true)
tic!(t,barrier=true)

# Parameters
ν=1.0
ρ=1.0
σ=1.0
B=VectorValue(0.0,0.0,10.0)
f=VectorValue(0.0,0.0,0.0)
u0=1.0
B0=norm(B)
L = 1.0 

# Reduced quantities
Re = u0*L/ν
Ha = B0*L*sqrt(σ/(ρ*ν))
N = Ha^2/Re
f̄ = (L/(ρ*u0^2))*f
B̄ = (1/B0)*B
α = 1.0
β = 1.0/Re
γ = N

# Domain and model
n=8
domain = (0,L,0,L,0,L)
cells = (n,n,n)
model = simplexify(CartesianDiscreteModel(domain,cells))
Ω = Interior(model)

# Boundary conditions
labels = get_face_labeling(model)
Γw = append!(collect(1:4),[9,10,13,14],collect(17:21),collect(23:26))
Γl = append!(collect(5:8),[11,12,15,16,22])
add_tag_from_tags!(labels,"wall",Γw)
add_tag_from_tags!(labels,"lid",Γl)
add_tag_from_tags!(labels,"insulating","boundary")
uw = VectorValue(0.0, 0.0, 0.0)
ul = VectorValue(1.0, 0.0, 0.0)
ji = VectorValue(0.0, 0.0, 0.0)

_params = Dict(
    :ptimer => t,
    :debug => false,
    :solve => true,
    :res_assemble => false,
    :jac_assemble => false,
    :model => model,
    :fluid => Dict(
      :domain=>model,
      :α=>α,
      :β=>β,
      :γ=>γ,
      :f=>f̄,
      :B=>B̄,
    ),
    :bcs => Dict(
        :u=>Dict(:tags=>["wall","lid"],:values=>[uw,ul]),
        :j=>Dict(:tags=>"insulating",:values=>ji)
    )
)
params = add_default_params(_params)
toc!(t,"pre_process")

tic!(t;barrier=true)
# ReferenceFEs
k = 2
T = Float64
model = params[:model]
D = num_cell_dims(model)
reffe_u = ReferenceFE(lagrangian,VectorValue{D,T},k)
reffe_p = ReferenceFE(lagrangian,T,k-1)
reffe_j = ReferenceFE(raviart_thomas,T,k-1)
reffe_φ = ReferenceFE(lagrangian,T,k-1)

# Test spaces
V_u = TestFESpace(model,reffe_u;dirichlet_tags=["wall","lid"])
V_p = TestFESpace(model,reffe_p;constraint=:zeromean)
V_j = TestFESpace(model,reffe_j;dirichlet_tags="insulating")
V_φ = TestFESpace(model,reffe_φ;conformity=:L2)
V = MultiFieldFESpace([V_u,V_p,V_j,V_φ])
#V   = MultiFieldFESpace([V_u,V_p,V_j,V_φ];style=BlockMultiFieldStyle())

# Trial spaces
U_u = TrialFESpace(V_u,[uw,ul])
U_j = TrialFESpace(V_j,ji)
U_p = TrialFESpace(V_p)
U_φ = TrialFESpace(V_φ)
U = MultiFieldFESpace([U_u,U_p,U_j,U_φ])
#U = MultiFieldFESpace([U_u,U_p,U_j,U_φ];style=BlockMultiFieldStyle())
toc!(t,"fe_spaces")

tic!(t;barrier=true)
# Weak form
#! ζ adds an Augmented-Lagragian term to both the preconditioner and teh weak form. 
#! Set to zero if not needed.
params[:ζ] = 0.0
res, jac = weak_form(params,k)
Tm = params[:matrix_type]
Tv = params[:vector_type]
assem = SparseMatrixAssembler(Tm,Tv,U,V)
op    = FEOperator(res,jac,U,V,assem)

xh = zero(U)
solver = NLSolver(show_trace=true,method=:newton)
xh,cache = solve!(xh,solver,op)
toc!(t,"solve")

tic!(t,barrier=true)
ūh,p̄h,j̄h,φ̄h = xh
uh = u0*ūh
ph = (ρ*u0^2)*p̄h
jh = (σ*u0*B0)*j̄h
φh = (u0*B0*L)*φ̄h
writevtk(Ω,"cavity",
  order=2,
  cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh])
toc!(t,"vtk")


info = Dict{Symbol,Any}()
info[:ncells] = num_cells(model)
info[:ndofs_u] = length(get_free_dof_values(ūh))
info[:ndofs_p] = length(get_free_dof_values(p̄h))
info[:ndofs_j] = length(get_free_dof_values(j̄h))
info[:ndofs_φ] = length(get_free_dof_values(φ̄h))
info[:ndofs] = length(get_free_dof_values(xh))
info[:Re] = Re
info[:Ha] = Ha
save("cavity.bson",info)


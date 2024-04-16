module TimeSteppingDev


using Gridap
using DrWatson
using GridapDistributed
using PartitionedArrays


using GridapMHD
using GridapMHD: default_solver_params
using GridapMHD: hunt_mesh
using GridapMHD: main
using GridapMHD: add_default_params

using Gridap.FESpaces
using Gridap.ODEs
using Gridap.ODEs.TransientFETools: TransientFEOperatorFromWeakForm
using Gridap.ODEs.TransientFETools: Nonlinear
using Gridap.ODEs.TransientFETools: rhs_error
import Gridap.ODEs.TransientFETools: TransientFEOperator

function TransientFEOperator(res::Function,jac::Function,jac_t::Function,
  trial,test,assem::Assembler)
  TransientFEOperatorFromWeakForm{Nonlinear}(res,rhs_error,(jac,jac_t),assem,(trial,∂t(trial)),test,1)
end

# for Δt in exp10.(-1:-1:-6)

Δt = 0.1
  println("Running Δt = $Δt")

  println("------------------------------------")

# Params
distribute = nothing
rank_partition=nothing
nc=(2,2,2)
ν=1.0
ρ=1.0
σ=1.0
B=(0.0,0.0,1.0)
# f=(0.0,0.0,1.0)
ζ=0.0
L=1.0
u0=1.0
B0=norm(VectorValue(B))
nsums = 10
vtk=true
title = "test"
path  = datadir()
debug = false
res_assemble = false
jac_assemble = false
solve = true
solver = :julia
verbose = true

t0 = 0.0
tf = Δt

info = Dict{Symbol,Any}()
params = Dict{Symbol,Any}(
  :debug=>debug,
  :solve=>solve,
  :res_assemble=>res_assemble,
  :jac_assemble=>jac_assemble,
)

# Manufactured solution [Zhang]

# CartesianDiscreteModel()
Re = κ = 1
B = VectorValue(0,0,1)
B0 = norm(B)

# # 2D
u(x,t::Real) = VectorValue(x[2]*exp(-t), x[1]*cos(t),0)
# p(x,t::Real) = sin(t)
# j(x,t::Real) = VectorValue(x[1]*sin(t),x[2]*cos(t),0)
# φ(x,t::Real) = 1

u(x,t::Real) = VectorValue(x[2]*(1-t), x[1]*(1-t),0)
p(x,t::Real) = 0
φ(x,t::Real) = 1

u(t::Real) = x -> u(x,t)
p(t::Real) = x -> p(x,t)
# j(t::Real) = x -> j(x,t)
φ(t::Real) = x -> φ(x,t)

#TODO: add f_j * v_j
j(x,t::Real) = u(x,t) × B - ∇(φ(t))(x)
j(t::Real) = x -> j(x,t)

f(x,t::Real) =
  (∂t(u))(x,t) +
  u(x,t) ⋅ ∇(u(t))(x) +
  - (1/Re) * Δ(u(t))(x) +
  ∇(p(t))(x) +
  - κ*( j(x,t) × B )

f(t::Real)  = x -> f(x,t)

# # 3D
# u(x,t) = VectorValue(x[3]*sin(t),x[1],x[2]*exp(-t))
# p(x,t) = 0
# j(x,t) = VectorValue(cos(t),t^2,0)
# φ(x,t) = 0

# Communicator
if isa(distribute,Nothing)
  @assert isa(rank_partition,Nothing)
  rank_partition = Tuple(fill(1,length(nc)))
  distribute = DebugArray
end
@assert length(rank_partition) == length(nc)
parts = distribute(LinearIndices((prod(rank_partition),)))

# Timer
t = PTimer(parts,verbose=verbose)
params[:ptimer] = t
tic!(t,barrier=true)

# Solver
if isa(solver,Symbol)
  solver = default_solver_params(Val(solver))
end
params[:solver] = solver

L/(ρ*u0^2)
# Reduced quantities
Re = u0*L/ν
Ha = B0*L*sqrt(σ/(ρ*ν))
N = Ha^2/Re
# f̄ = t -> x -> (L/(ρ*u0^2))*f(t)(x)
# B̄ = (1/B0)*B
f̄ = f
B̄ = B
α = 1.0
β = 1.0/Re
γ = N

# DiscreteModel in terms of reduced quantities

nc
np = rank_partition
domain = ( 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 )
model = CartesianDiscreteModel(parts,np,domain,nc)
params[:model] = model

using Gridap.ReferenceFEs
using Gridap.Geometry
D = 3
d = 2
pt = HEX
D = num_dims(pt)
facets = 2*(D-d)+1:2*(D)
faces = facets .+ get_offset(pt,D-1)
ids = get_faces(pt)[faces]
ids = reduce(union,ids) |> sort
ids = collect(faces)

map(local_views(model)) do model
  labels = get_face_labeling(model)
  add_tag!(labels,"boundary_2D",ids)
end

facets = 1:2
faces = facets .+ get_offset(pt,D-1)
ids = collect(faces)
map(local_views(model)) do model
  labels = get_face_labeling(model)
  add_tag!(labels,"top_bottom",ids)
end


writevtk(model,"model")

Ω = Interior(model)
if debug && vtk
  writevtk(model,"data/manufactured_model")
end

map(num_dims,local_views(model))

params[:fluid] = Dict(
  :domain=>nothing,
  :α=>α,
  :β=>β,
  :γ=>γ,
  :f=>f̄,
  :B=>B̄,
  :ζ=>ζ,
)

# Boundary conditions

params[:bcs] = Dict(
  :u=>Dict(:tags=>"boundary",:values=>u),
  :j=>Dict(:tags=>"boundary_2D",:values=>j),
  :φ=>Dict(:domain=>"top_bottom",:value=>φ),
)

params[:fespaces] = Dict(
  :k => 2,
  :p_constraint => :zeromean)

toc!(t,"pre_process")

ode_solver_params = Dict(:θ=>0.5)

initial_values = Dict(
  :u => u(0),
  :p => p(0),
  :j => j(0),
  :φ => φ(0),
)

params[:ode] = Dict(
  :solver => :theta,
  :t0 => t0,
  :tf => tf,
  :Δt => Δt,
  :solver_params => ode_solver_params,
  :U0 => :value,
  :initial_values => initial_values,
)

xh,fullparams,info = main(params;output=info)


# Bug implementing Dirichlet BC

# Post process

if L == 1.0
  Ω_phys = Ω
else
  Ω_phys = _warp(model,Ω,L)
end

results = Dict(
  :Δt => Δt,
  :t0 => t0,
  :tf => tf,
  :t => Δt:Δt:tf,
  :el2_uh => [],
  :eh1_uh => [],
  :el2_jh => [],
  :el2_ph => [],
  :el2_φh => [],
)

tic!(t,barrier=true)
pvd = createpvd(parts,"transient")

dΩ = Measure(Ω_phys,2)
ℓ2(u,dΩ) = √(∑( ∫(u ⊙ u)dΩ ))
h1(u,dΩ) = √(∑( ∫( ∇(u) ⊙ ∇(u) + u ⊙ u)dΩ ))
t0 = 0

pvd[t0] = createvtk(Ω_phys,"transient_0",
  order=2,
  cellfields=[
    "uh"=>u(t0),"ph"=>p(t0),"jh"=>j(t0),"phi"=>φ(t0),
    "u"=>u(t0),"j"=>j(t0)])

for (i,(xht,t)) in enumerate(xh)
  ūh,p̄h,j̄h,φ̄h = xht
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (u0*B0*L)*φ̄h

  e_uh = uh - u(t)
  e_jh = jh - j(t)
  e_ph = (ph) - (p(t))
  e_φh = (φh) - (φ(t))

  e_uh_l2 = ℓ2(e_uh,dΩ)
  e_jh_l2 = ℓ2(e_jh,dΩ)
  e_ph_l2 = ℓ2(e_ph,dΩ)
  e_φh_l2 = ℓ2(e_φh,dΩ)

  e_uh_h1 = h1(e_uh,dΩ)

  @show e_uh_l2
  @show e_uh_h1
  @show e_jh_l2
  @show e_ph_l2
  @show e_φh_l2

  push!(results[:el2_uh],e_uh_l2)
  push!(results[:eh1_uh],e_uh_h1)
  push!(results[:el2_jh],e_jh_l2)
  push!(results[:el2_ph],e_ph_l2)
  push!(results[:el2_φh],e_φh_l2)

  pvd[t] = createvtk(Ω_phys,"transient_$i",
    order=2,
    cellfields=[
    "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
    "u"=>u(t),"j"=>j(t)])
end
savepvd(pvd)
toc!(t,"time_stepping")

save("results_$Δt.jld2",tostringdict(results))

# TODO:

# Convergence test
# end

end # module

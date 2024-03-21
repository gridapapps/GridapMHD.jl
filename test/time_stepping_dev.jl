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

# Params
distribute = nothing
rank_partition=nothing
nc=(4,4)
ν=1.0
ρ=1.0
σ=1.0
B=(0.0,10.0,0.0)
f=(0.0,0.0,1.0)
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
BL_adapted = true
kmap_x = 1
kmap_y = 1
ranks_per_level = nothing

nc=(10,10)
Δt = 0.2


info = Dict{Symbol,Any}()
params = Dict{Symbol,Any}(
  :debug=>debug,
  :solve=>solve,
  :res_assemble=>res_assemble,
  :jac_assemble=>jac_assemble,
)

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

# Reduced quantities

Re = u0*L/ν
Ha = B0*L*sqrt(σ/(ρ*ν))
N = Ha^2/Re
f̄ = (L/(ρ*u0^2))*VectorValue(f)
B̄ = (1/B0)*VectorValue(B)
f̄ₜ = t -> f̄*t
B̄ₜ = t -> B̄*t
α = 1.0
β = 1.0/Re
γ = N

# DiscreteModel in terms of reduced quantities

model = hunt_mesh(parts,params,nc,rank_partition,L,Ha,kmap_x,kmap_y,BL_adapted,ranks_per_level)
Ω = Interior(model)
if debug && vtk
  writevtk(model,"data/hunt_model")
end

params[:fluid] = Dict(
  :domain=>nothing,
  :α=>α,
  :β=>β,
  :γ=>γ,
  :f=>f̄ₜ,
  :B=>B̄ₜ,
  :ζ=>ζ,
)

# Boundary conditions

params[:bcs] = Dict(
  :u=>Dict(:tags=>"noslip"),
  :j=>Dict(:tags=>"insulating"),
)

toc!(t,"pre_process")

ode_solver_params = Dict(:θ=>0.5)

params[:ode] = Dict(
  :solver => :theta,
  :t0 => 0.0,
  :tf => 1.0,
  :Δt => Δt,
  :solver_params => ode_solver_params,
  :X0 => :zero,
)

xh,fullparams,info = main(params;output=info)

  # Post process

  using GridapMHD: analytical_hunt_u, analytical_hunt_j

  μ = ρ*ν
  grad_pz = -f[3]/ρ
  u(x) = analytical_hunt_u(L,L,μ,grad_pz,Ha,nsums,x)
  j(x) = analytical_hunt_j(L,L,σ,μ,grad_pz,Ha,nsums,x)
  u_ref(x) = analytical_hunt_u(L,L,μ,grad_pz,Ha,2*nsums,x)
  j_ref(x) = analytical_hunt_j(L,L,σ,μ,grad_pz,Ha,2*nsums,x)

if L == 1.0
  Ω_phys = Ω
else
  Ω_phys = _warp(model,Ω,L)
end

tic!(t,barrier=true)
pvd = createpvd(parts,"hunt_transient")
for (i,(xht,t)) in enumerate(xh)
  ūh,p̄h,j̄h,φ̄h = xht
  uh = u0*ūh
  ph = (ρ*u0^2)*p̄h
  jh = (σ*u0*B0)*j̄h
  φh = (u0*B0*L)*φ̄h
  pvd[t] = createvtk(Ω_phys,"hunt_transient_$i",
    # order=2,
    cellfields=[
    "uh"=>uh,"ph"=>ph,"jh"=>jh,"phi"=>φh,
    "u"=>u,"j"=>j,"u_ref"=>u_ref,"j_ref"=>j_ref])
end
savepvd(pvd)
toc!(t,"time_stepping")


# TODO:
#
# Setup transient benchmark
# Setup initial values
# Plot initial values
# Convergence test


end # module

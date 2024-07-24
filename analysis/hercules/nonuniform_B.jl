using DrWatson
include("scripts/mpi_scripts.jl")

setup = "
using SparseMatricesCSR
using GridapPETSc
using GridapMHD: snes_postpro"
solver = "Dict(
  :solver => :petsc,
  :matrix_type    => SparseMatrixCSR{0,PetscScalar,PetscInt},
  :vector_type    => Vector{PetscScalar},
  :solver_postpro => ((cache,info) -> snes_postpro(cache,info)),
  :petsc_options  => \"-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_28 1 -mat_mumps_icntl_29 2 -mat_mumps_icntl_4 3 -mat_mumps_cntl_1 0.001\",
  :niter          => 100,
  :rtol           => 1e-5,
)"

function run_channel(dict)
  run_channel(;dict2ntuple(dict)...)
end

function run_channel(;force=false,L=1,w=0.1,kwargs...)

  @unpack np,nc,Ha,Re,niter,convection,nonuniform_B,γ = kwargs

  N = prod(nc)
  prefix = "channel_N$N"

  jobname = savename(prefix,kwargs,equals="",ignores=["np"])
  pvtu = datadir("$jobname.pvtu")
  if !force && isfile(pvtu)
    @info "Skipping $jobname"
    return
  else
    @info "Running $jobname"
  end

  run_driver(
    app="channel",
    setup=setup,
    nc=nc,
    np=np,
    backend=:mpi,
    sizes=(L,w,w),
    ν=1,
    B0=Ha/w,
    u0=Re/w,
    bl_orders=(1,2,2),
    solver=":petsc",
    niter=10,
    inlet=:shercliff,
    initial_value=:zero,
    convection=convection,
    nonuniform_B=nonuniform_B,
    γB=γ,
    title=jobname,
    vtk=true)
end


params = Dict(
  :Ha => [100],
  :Re => [10],
  :niter => nothing,
  :np => (4,2,2),
  :nc => [(40,20,10)],
  L=30,
  w=2,
  :nonuniform_B => true,
  :γ => 0.45,
  :convection => true,
  :initial_value => :inlet,
  :force => true,
)

for p in dict_list(params)
  run_channel(p)
end

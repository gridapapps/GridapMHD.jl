
using DrWatson
include("scripts/mpi_scripts.jl")

setup = "
using GridapMHD: expansion
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

function run_expansion(dict)
  run_expansion(;dict2ntuple(dict)...)
end

function run_expansion(;force=false,kwargs...)

  prefix = "expansion"

  jobname = savename(prefix,kwargs,equals="",ignores=["np"])
  @unpack np,mesh,Ha,Re,niter,convection = kwargs
  pvtu = datadir("$jobname.pvtu")
  if !force && isfile(pvtu)
    @info "Skipping $jobname"
    return
  else
    @info "Running $jobname"
  end

  run_driver(
    app="expansion",
    setup=setup,
    mesh=mesh,
    np=np,
    backend=:mpi,
    Ha = Ha,
    N = Ha^2/Re,
    cw = 0.0,
    Z = 4.0,
    b = 1.0,
    inlet=:shercliff,
    debug=false,
    vtk=true,
    title=jobname,
    solver=solver,
    initial_value=:zero,
    niter=niter,
    convection=convection)
end



params = Dict(
  :Ha => [10,100],
  :Re => 1,
  :niter => 1,
  :np => 16,
  :mesh => ["710","6k"],
  :convection => false,
)


for p in dict_list(params)
  run_expansion(p)
end

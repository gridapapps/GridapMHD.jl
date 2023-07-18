using GridapMHD: hunt
using GridapMHD: expansion

hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  nsums=50,
  path=".",
  nruns=1,
  verbose=true,
  title="warmup_hunt_julia",
  solver=:julia,
  BL_adapted = true 
)

hunt(
  nc=(4,4),
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  nsums=50,
  path=".",
  nruns=1,
  verbose=true,
  title="warmup_hunt_petsc",
  solver=:petsc,
  BL_adapted = true,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_4 2  -mat_mumps_icntl_28 2  -mat_mumps_icntl_29 2"
)

hunt(
  nc=(4,4),
  np=(1,1),
  backend=:mpi,
  L=1.0,
  B=(0.,50.,0.),
  debug=false,
  vtk=true,
  nsums=50,
  path=".",
  nruns=1,
  verbose=true,
  title="warmup_hunt_mpi",
  solver=:petsc,
  BL_adapted = true,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_4 2  -mat_mumps_icntl_28 2  -mat_mumps_icntl_29 2"
)

expansion(
  mesh="710",
  Ha = 10.0,
  N = 10.0,
  cw = 0.01,
  debug=false,
  vtk=true,
  title="warmup_gmsh",
  solver=:julia
)

expansion(
  mesh="710",
  Ha = 10.0,
  N = 10.0,
  cw = 0.01,
  debug=false,
  vtk=true,
  title="warmup_gmsh_petsc",
  solver=:petsc,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_4 2  -mat_mumps_icntl_28 2  -mat_mumps_icntl_29 2"
)


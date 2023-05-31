using GridapMHD: hunt

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
  title="warmup_hunt",
  solver=:petsc,
  kmap=1,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
  res_assemble=true,
  jac_assemble=true,
  solve=true
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
  kmap=1,
  petsc_options="-snes_monitor -ksp_error_if_not_converged true -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps",
  res_assemble=true,
  jac_assemble=true,
  solve=true
)

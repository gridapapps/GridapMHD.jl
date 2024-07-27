using MPI

if !MPI.Initialized()
  MPI.Init()
end

using GridapMHD: hunt
using GridapMHD: expansion
using GridapMHD: pipe

hunt(
  title="warmup_hunt_res_jac",
  solve=false,
  res_assemble=true,
  jac_assemble=true
)

hunt(
  title="warmup_hunt_julia",
  solver=:julia
)

hunt(
  title="warmup_hunt_petsc",
  solver=:petsc
)

hunt(
  np=(2,2),
  backend=:sequential,
  title="warmup_hunt_np4",
  solver=:petsc
)

hunt(
  np=(1,1),
  backend=:mpi,
  title="warmup_hunt_mpi",
  solver=:julia
)

expansion(title="warmup_gmsh")
expansion(inlet=:shercliff,title="warmup_gmsh_shercliff")
expansion(np=2,backend=:sequential,title="warmup_gmsh_np2")
expansion(np=1,backend=:mpi,title="warmup_gmsh_mpi")

hunt(
  title="warmup_pipe_julia",
  solver=:julia
)

hunt(
  title="warmup_pipe_petsc",
  solver=:petsc
)

hunt(
  np=(1,1),
  backend=:mpi,
  title="warmup_pipe_mpi",
  solver=:julia
)

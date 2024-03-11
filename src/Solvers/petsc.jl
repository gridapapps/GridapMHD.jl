
get_block_solver(::Val{:petsc_mumps},params)         = PETScLinearSolver(petsc_mumps_setup)
get_block_solver(::Val{:petsc_amg},params)           = PETScLinearSolver(petsc_amg_setup)
get_block_solver(::Val{:petsc_cg_jacobi},params)     = PETScLinearSolver(petsc_cg_setup)
get_block_solver(::Val{:petsc_gmres_schwarz},params) = PETScLinearSolver(petsc_gmres_schwarz_setup)
get_block_solver(::Val{:petsc_gmres_amg},params)     = PETScLinearSolver(petsc_gmres_amg_setup)
get_block_solver(::Val{:petsc_from_options},params)  = PETScLinearSolver()

function petsc_mumps_setup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  #@check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  7, 0)
  #@check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  #@check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
  #@check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-6)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

# Standalone AMG solver, 5 iterations
function petsc_amg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(5)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPRICHARDSON)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

# GMRES + Additive Schwartz preconditioner
function petsc_gmres_schwarz_setup(ksp)
  rtol = PetscScalar(1.e-6)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  # GMRES solver
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)

  # Additive Schwartz preconditioner
  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCASM)

  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

# CG + Jacobi preconditioner, 10 iterations
function petsc_cg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(10)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPCG)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCJACOBI)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

function petsc_gmres_amg_setup(ksp)
  rtol = PetscScalar(1.e-6)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

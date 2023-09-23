
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using Gridap.ReferenceFEs

function mumps_setup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
  # percentage increase in the estimated working space
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-6)
end

function gmres_schwarz_setup(ksp)
  rtol = PetscScalar(1.e-5)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  # GMRES solver
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)

  # Additive Schwartz preconditioner
  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCASM)
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

np = (2,2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

# Geometry
n = 20
model = CartesianDiscreteModel(ranks,np,(0,1,0,1,0,1),(n,n,3))

labels = get_face_labeling(model);
tags_j = append!(collect(1:20),[25,26])
add_tag_from_tags!(labels,"insulating",tags_j)
j_bc(x) = VectorValue(0.0,0.0,0.0)

# FESpaces
k = 2
reffe_j = ReferenceFE(raviart_thomas,Float64,k-1)
V_j = TestFESpace(model,reffe_j;dirichlet_tags="insulating")
U_j = TrialFESpace(V_j,j_bc)

γ = 1.0
Ω = Interior(model)
dΩ = Measure(Ω,2*k)
a_j(j,v_j) = ∫(γ*j⋅v_j + γ*(∇⋅j)⋅(∇⋅v_j))*dΩ
D_j = assemble_matrix(a_j,U_j,V_j)

options = """
  -ksp_type gmres
  -ksp_rtol 1.0e-5
  -ksp_atol 1.0e-14
  -ksp_monitor
  -ksp_converged_reason
  -pc_type asm
  -pc_asm_overlap 10
  -pc_asm_type restrict
  -pc_asm_blocks 64
  -sub_ksp_type preonly
  -sub_pc_type lu
"""

GridapPETSc.with(args=split(options)) do 
  solver = PETScLinearSolver()
  test_solver(solver,D_j)
end

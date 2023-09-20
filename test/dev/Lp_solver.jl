
using PartitionedArrays
using Gridap, GridapPETSc, GridapSolvers, GridapDistributed

using GridapSolvers.LinearSolvers
using Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.Geometry

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
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  14, 1000)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[], 3, 1.0e-6)
end

function gmres_amg_setup(ksp)
  rtol = PetscScalar(1.e-5)
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = GridapPETSc.PETSC.PETSC_DEFAULT

  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
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

function amg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(10)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPRICHARDSON)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  #@check_error_code GridapPETSc.PETSC.PCGAMGSetType(pc[],GridapPETSc.PETSC.PCGAMGAGG) # or PCGAMGCLASSICAL
  #@check_error_code GridapPETSc.PETSC.PCGAMGSetNSmooths(pc[],0)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
end

get_edge_measures(Ω::Triangulation,dΩ) = CellField(get_array(∫(1)dΩ),Ω)
function get_edge_measures(Ω::GridapDistributed.DistributedTriangulation,dΩ) 
  return CellField(map(get_array,local_views(∫(1)*dΩ)),Ω)
end

np = (2,2,1)
ranks = with_debug() do distribute
  distribute(LinearIndices((prod(np),)))
end

# Geometry
n = 16
model = CartesianDiscreteModel(ranks,np,(0,1,0,1,0,1),(n,n,n))

# FESpaces
k = 2
p_ref(x) = 0.0
reffe_p = ReferenceFE(lagrangian,Float64,k-1,space=:P)
V_p = TestFESpace(model,reffe_p)#;dirichlet_tags=["boundary"])
U_p = TrialFESpace(V_p,p_ref)

Ω = Interior(model)
Γ = Boundary(model)
Λ = Skeleton(model)

dΩ = Measure(Ω,2*k)
dΓ = Measure(Γ,2*k)
dΛ = Measure(Λ,2*k)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

h_e_Λ = get_edge_measures(Λ,dΛ)
h_e_Γ = get_edge_measures(Γ,dΓ)

a_cg(p,v_p) = ∫(∇(p)⋅∇(v_p))*dΩ

β = 50.0
aΛ(u,v) = ∫(-jump(u⋅n_Λ)⋅mean(∇(v)) - mean(∇(u))⋅jump(v⋅n_Λ))*dΛ + ∫(β/h_e_Λ*jump(u⋅n_Λ)⋅jump(v⋅n_Λ))*dΛ
aΓ(u,v) = ∫(-(∇(u)⋅n_Γ)⋅v - u⋅(∇(v)⋅n_Γ))*dΓ + ∫(β/h_e_Γ*(u⋅n_Γ)⋅(v⋅n_Γ))*dΓ
a_dg(p,v_p) = a_cg(p,v_p) + aΛ(p,v_p) + aΓ(p,v_p)

Δp = assemble_matrix(a_dg,U_p,V_p)

# Solver
b = GridapSolvers.allocate_col_vector(Δp)
x = GridapSolvers.allocate_col_vector(Δp)

fill!(x,0.0)
fill!(b,1.0)
@time begin
  GridapPETSc.with(args=split("-ksp_converged_reason")) do 
    solver = PETScLinearSolver(gmres_schwarz_setup)
    ns = numerical_setup(symbolic_setup(solver, Δp), Δp)
    solve!(x,ns,b)
    println(norm(b - Δp*x))
  end
end


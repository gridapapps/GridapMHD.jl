
using GridapPETSc: PETScMatrix, PETScVector, PETScLinearSolverSS
using GridapPETSc.PETSC
using MPI

mutable struct CachedPETScLinearSolverNS <: NumericalSetup
  B::PETScMatrix
  ksp::Ref{KSP}
  initialized::Bool
  caches
  function CachedPETScLinearSolverNS(B::PETScMatrix,caches)
    new(B,Ref{KSP}(),false,caches)
  end
end

function Init(a::CachedPETScLinearSolverNS)
  @assert Threads.threadid() == 1
  GridapPETSc._NREFS[] += 1
  a.initialized = true
  finalizer(Finalize,a)
end

function Finalize(ns::CachedPETScLinearSolverNS)
  if ns.initialized && GridapPETSc.Initialized()
    if ns.B.comm == MPI.COMM_SELF
      @check_error_code PETSC.KSPDestroy(ns.ksp)
    else
      @check_error_code PETSC.PetscObjectRegisterDestroy(ns.ksp[].ptr)
    end
    ns.initialized = false
    @assert Threads.threadid() == 1
    GridapPETSc._NREFS[] -= 1
  end
  nothing
end

function get_petsc_caches(A::AbstractMatrix)
  X = convert(PETScVector,allocate_col_vector(A))
  B = convert(PETScVector,allocate_col_vector(A))
  return X,B
end

function Algebra.numerical_setup(ss::PETScLinearSolverSS,A::AbstractMatrix)
  B = convert(PETScMatrix,A)
  caches = get_petsc_caches(A)
  ns = CachedPETScLinearSolverNS(B,caches)
  @check_error_code PETSC.KSPCreate(B.comm,ns.ksp)
  @check_error_code PETSC.KSPSetOperators(ns.ksp[],ns.B.mat[],ns.B.mat[])
  ss.solver.setup(ns.ksp)
  @check_error_code PETSC.KSPSetUp(ns.ksp[])
  Init(ns)
end

function Algebra.solve!(x::AbstractVector{PetscScalar},ns::CachedPETScLinearSolverNS,b::AbstractVector{PetscScalar})
  X, B = ns.caches
  copy!(B,b)
  @check_error_code PETSC.KSPSolve(ns.ksp[],B.vec[],X.vec[])
  copy!(x,X)
  return x
end

function Algebra.numerical_setup!(ns::CachedPETScLinearSolverNS,A::AbstractMatrix)
  copy!(ns.B,A)
  @check_error_code PETSC.KSPSetOperators(ns.ksp[],ns.B.mat[],ns.B.mat[])
  @check_error_code PETSC.KSPSetUp(ns.ksp[])
  return ns
end

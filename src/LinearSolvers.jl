module LinearSolvers

using Gridap
using Gridap.Algebra
using IterativeSolvers
using IncompleteLU

import Gridap.Algebra: LinearSolver
import Gridap.Algebra: symbolic_setup
import Gridap.Algebra: numerical_setup
import Gridap.Algebra: numerical_setup!
import Gridap.Algebra: solve!

export GmresSolver
export symbolic_setup
export numerical_setup
export numerical_setup!
export solve!


struct GmresSolver <: LinearSolver
  initially_zero::Bool
  preconditioner
  precond_kwargs::Dict
  function GmresSolver(;
    initially_zero=false,preconditioner=IterativeSolvers.Identity,
    precond_kwargs...)
    new(initially_zero,preconditioner,precond_kwargs)
  end
end

struct GmresSymbolicSetup <: SymbolicSetup
  initially_zero::Bool
  preconditioner
  precond_kwargs::Dict
end

mutable struct GmresNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
  A::T
  initially_zero::Bool
  preconditioner
  precond_kwargs::Dict
end

function symbolic_setup(s::GmresSolver, mat::AbstractMatrix)
  GmresSymbolicSetup(s.initially_zero, s.preconditioner, s.precond_kwargs)
end

function numerical_setup(ss::GmresSymbolicSetup,mat::AbstractMatrix)
  GmresNumericalSetup(mat, ss.initially_zero,
    ss.preconditioner, ss.precond_kwargs)
end

function numerical_setup!(ns::GmresNumericalSetup, mat::AbstractMatrix)
  ns.A = mat
end

function get_preconditioner(ns::GmresNumericalSetup)

end

function solve!(
  x::AbstractVector,ns::GmresNumericalSetup,b::AbstractVector)
  p=ns.preconditioner(ns.A; ns.precond_kwargs...)
  # ilu(ns.A, τ=ns.τ)
  gmres!(x, ns.A, b,verbose=true, Pl=p, initially_zero=ns.initially_zero)
end

end # module

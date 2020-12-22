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
  preconditioner
  precond_kwargs::Dict
  verbose::Bool
  function GmresSolver(;
    verbose=false, preconditioner=IterativeSolvers.Identity, precond_kwargs...)
    new(preconditioner,precond_kwargs, verbose)
  end
end

struct GmresSymbolicSetup <: SymbolicSetup
  preconditioner
  precond_kwargs::Dict
  verbose::Bool
end

mutable struct GmresNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
  A::T
  preconditioner
  precond_kwargs::Dict
  verbose::Bool
end

function symbolic_setup(s::GmresSolver, mat::AbstractMatrix)
  GmresSymbolicSetup(s.preconditioner, s.precond_kwargs, s.verbose)
end

function numerical_setup(ss::GmresSymbolicSetup,mat::AbstractMatrix)
  GmresNumericalSetup(mat, ss.preconditioner, ss.precond_kwargs, ss.verbose)
end

function numerical_setup!(ns::GmresNumericalSetup, mat::AbstractMatrix)
  ns.A = mat
end

function get_preconditioner(ns::GmresNumericalSetup)

end

function solve!(
  x::AbstractVector,ns::GmresNumericalSetup,b::AbstractVector)
  p=ns.preconditioner(ns.A; ns.precond_kwargs...)
  # initialize the solution to 0
  x .= 0.0
  gmres!(x, ns.A, b,verbose=ns.verbose, Pl=p, initially_zero=true)
end

end # module

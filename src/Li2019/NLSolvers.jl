
function Gridap.Algebra.solve!(x::AbstractVector,nls::NewtonRaphsonSolver,op::NonlinearOperator,cache::Nothing)
  b  = residual(op, x)
  A  = jacobian(op, x)
  dx = allocate_col_vector(A)
  ns = numerical_setup(symbolic_setup(nls.ls, A), A)

  Gridap.Algebra._solve_nr!(x,A,b,dx,ns,nls,op)
  return Gridap.Algebra.NewtonRaphsonCache(A,b,dx,ns)
end

function Gridap.Algebra._solve_nr!(x,A,b,dx,ns,nls,op)
  
  # Check for convergence on the initial residual
  isconv, conv0 = Gridap.Algebra._check_convergence(nls,b)
  if isconv; return; end

  # Newton-like iterations
  for nliter in 1:nls.max_nliters

    # Solve linearized problem
    rmul!(b,-1)
    solve!(dx,ns,b)
    x .+= dx

    # Check convergence for the current residual
    residual!(b, op, x)
    isconv = Gridap.Algebra._check_convergence(nls, b, conv0)
    if isconv; return; end

    if nliter == nls.max_nliters
      @unreachable
    end

    # Assemble jacobian (fast in-place version)
    # and prepare solver
    jacobian!(A, op, x)
    numerical_setup!(ns,A,x)
  end
end

function Gridap.Algebra._check_convergence(nls,b)
  parts = GridapDistributed.get_parts(b)
  m0 = maximum(abs,b)
  i_am_main(parts) && println(" >> Non-linear residual: e_a = $m0, e_r = 1.0")
  return (false, m0)
end

function Gridap.Algebra._check_convergence(nls,b,m0)
  parts = GridapDistributed.get_parts(b)
  m = maximum(abs,b)
  i_am_main(parts) && println(" >> Non-linear residual: e_a = $m, e_r = $(m/m0)")
  return m < nls.tol * m0
end

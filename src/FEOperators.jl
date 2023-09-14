"""

"""
struct FEOperatorMHD <: Gridap.FESpaces.FEOperator
  op_global :: FEOperator # global operator
  op_u      :: FEOperator # operator for u
end

Gridap.FESpaces.get_test(op::FEOperatorMHD)  = get_test(op.op_global)
Gridap.FESpaces.get_trial(op::FEOperatorMHD) = get_trial(op.op_global)

# Delegate most of the methods to the global operator
Gridap.FESpaces.allocate_residual(op::FEOperatorMHD,uh) = allocate_residual(op.op_global,uh)
Gridap.FESpaces.residual!(b::AbstractVector,op::FEOperatorMHD,uh) = residual!(b,op.op_global,uh)
Gridap.FESpaces.residual_and_jacobian(op::FEOperatorMHD,uh) = residual_and_jacobian(op.op_global,uh)
Gridap.FESpaces.residual_and_jacobian!(b::AbstractVector,A::AbstractMatrix,op::FEOperatorMHD,uh) = residual_and_jacobian!(b,A,op.op_global,uh)

# When updating the matrix, only update the u-u block
Gridap.FESpaces.allocate_jacobian(op::FEOperatorMHD,uh) = allocate_jacobian(op.op_global,uh)
Gridap.FESpaces.jacobian(op::FEOperatorMHD,uh) = jacobian(op.op_global,uh)
function Gridap.FESpaces.jacobian!(A::AbstractBlockMatrix,op::FEOperatorMHD,xh)
  uh, ph, jh, φh = xh
  jacobian!(A[Block(1,1)],op.op_u,uh)
  return A
end

function weakform_uu(params,k)
  fluid = params[:fluid]
  Ωf  = _interior(params[:model],fluid[:domain])
  dΩf = Measure(Ωf,2*k)

  α  = fluid[:α]
  β  = fluid[:β]
  f  = fluid[:f]
  ζ  = params[:ζ]

  if ζ !== nothing 
    a(u,dv) = ∫(β*(∇(u)⊙∇(dv)) ) * dΩf
  else
    a(u,dv) = ∫(β*(∇(u)⊙∇(dv)) + ζ*(∇⋅u)*(∇⋅v)) * dΩf
  end
  l(dv) = ∫( dv⋅f )*dΩf
  c(u,dv) = ∫( α*dv⋅(conv∘(u,∇(u))) ) * dΩf
  dc(u,du,dv) = ∫( α*dv⋅( (conv∘(u,∇(du))) + (conv∘(du,∇(u))) ) ) * dΩf

  res(u,dv) = c(u,dv) + a(u,dv) - l(dv)
  jac(u,du,dv) = dc(u,du,dv) + a(du,dv)
  return res, jac
end

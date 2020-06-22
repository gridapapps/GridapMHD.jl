module InductionlessMHD

using Gridap

function dimensionless_a(X,Y,Re,N,B)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  (1/Re)*inner(∇(u),∇(v_u)) - p*(∇⋅v_u) - N*(j×B)⋅v_u +
  (∇⋅u)*v_p +
  j⋅v_j - φ*(∇⋅v_j) - (u×B)⋅v_j +
  (∇⋅j)*v_φ
end

function dimensionless_l(Y,f_u)
  v_u, v_p, v_j, v_φ = Y

  v_u⋅f_u
end

@law conv(u,∇u) = (∇u')⋅u
@law dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

function convective_term(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  v_u⋅conv(u,∇(u))
end

function dconvective_term(X,dX,Y)
  u  , p  , j  , φ   = X
  du , dp , dj , dφ  = dX
  v_u, v_p, v_j, v_φ = Y

  v_u⋅dconv(du,∇(du),u,∇(u))
end

function dimensionless_residual(X,Y,Re,N,B)
  dimensionless_a(X,Y,Re,N,B) +
  convective_term(X,Y)
end

function dimensionless_residual(X,Y,Re,N,B,f_u)
  dimensionless_a(X,Y,Re,N,B) +
  convective_term(X,Y) -
  dimensionless_l(Y,f_u)
end

function dimensionless_jacobian(X,dX,Y,Re,N,B)
  dimensionless_a(dX,Y,Re,N,B) +
  dconvective_term(X,dX,Y)
end

function dimensionless_darcy_l_Γ(Y,n,g_φ)
  v_u, v_p, v_j, v_φ = Y

  -(v_j⋅n)*g_φ
end

function dimensionless_conducting_wall(X,Y,n,c_w)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y
  # (v_j⋅j)*0.0
  (v_j⋅n)*(j⋅n) - c_w * ((v_j⋅n)*(n⋅∇(j)⋅n))
end

end # module

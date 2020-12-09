module InductionlessMHD

using Gridap

function dimensionless_a(X,Y,Re,N,B,dΩ)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  ∫( (1/Re)*inner(∇(u),∇(v_u)) - p*(∇⋅v_u) - N*(j×B)⋅v_u +
     (∇⋅u)*v_p +
     j⋅v_j - φ*(∇⋅v_j) - (u×B)⋅v_j +
     (∇⋅j)*v_φ ) * dΩ
end

function dimensionless_l(Y,f_u,dΩ)
  v_u, v_p, v_j, v_φ = Y

  ∫( v_u⋅f_u ) * dΩ
end

conv(u,∇u) = (∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

function convective_term(X,Y,dΩ)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  ∫( v_u⋅conv(u,∇(u)) ) * dΩ
end

function dconvective_term(X,dX,Y,dΩ)
  u  , p  , j  , φ   = X
  du , dp , dj , dφ  = dX
  v_u, v_p, v_j, v_φ = Y

  ∫( v_u⋅dconv(du,∇(du),u,∇(u)) ) * dΩ
end

function dimensionless_residual(X,Y,Re,N,B,dΩ)
  dimensionless_a(X,Y,Re,N,B,dΩ) +
  convective_term(X,Y,dΩ)
end

function dimensionless_residual(X,Y,Re,N,B,f_u,dΩ)
  dimensionless_a(X,Y,Re,N,B,dΩ) +
  convective_term(X,Y,dΩ) -
  dimensionless_l(Y,f_u,dΩ)
end

function dimensionless_jacobian(X,dX,Y,Re,N,B,dΩ)
  dimensionless_a(dX,Y,Re,N,B,dΩ) +
  dconvective_term(X,dX,Y,dΩ)
end

function dimensionless_darcy_l_Γ(Y,n,g_φ,dΓ)
  v_u, v_p, v_j, v_φ = Y

  ∫( -(v_j⋅n)*g_φ ) * dΓ
end

function dimensionless_conducting_wall(X,Y,n,c_w,dΓ;α=1.0)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y
  ∫( α*((v_j⋅n)*(j⋅n) + c_w * (v_j⋅n)*(n⋅(∇(j)⋅n))) ) * dΓ
end

end # module

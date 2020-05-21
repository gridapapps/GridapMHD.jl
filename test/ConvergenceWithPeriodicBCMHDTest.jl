module ConvergenceWithPeriodicBCMHDTest

include("../src/GridapMHD.jl")
using .GridapMHD

include("../src/Defaults.jl")

using Gridap
using Test
using Polynomials: fit


u(x) = VectorValue(one(x[1]),x[3]^4+x[2]^4,-x[3]^4-x[2]^4)
p(x) = sin(x[3]+x[2])
j(x) = VectorValue(one(x[1]),(x[3])+x[2]^4,-x[3]^4-x[2])
φ(x) = sin(x[3]+x[2])
B = VectorValue(1.0,1.0,1.0)

∇u(x) = ∇(u)(x)
Δu(x) = Δ(u)(x)
∇p(x) = ∇(p)(x)
∇φ(x) = ∇(φ)(x)

function ∇u_n(x)
  if (x[2] < -0.49999)
    n = VectorValue(0.0,-1.0,0.0)
  elseif (x[2] > 0.49999)
    n = VectorValue(0.0,1.0,0.0)
  elseif (x[3] < -0.49999)
    n = VectorValue(0.0,0.0,-1.0)
  elseif (x[3] > 0.49999)
    n = VectorValue(0.0,0.0,1.0)
  elseif (x[1] < -0.49999)
    n = VectorValue(-1.0,0.0,0.0)
  else
    n = VectorValue(1.0,0.0,0.0)
  end
  n*∇u(x)
end

f_u(x) = (∇u(x)')*u(x) - Δu(x) + ∇p(x) - vprod(j(x),B)
f_p(x) = (∇*u)(x)
f_j(x) = j(x) + ∇φ(x) - vprod(u(x),B)
f_φ(x) = (∇*j)(x)

g_u(x) = u(x)
g_j(x) = j(x)

eu_l2 = Vector{Float64}()
eu_h1 = Vector{Float64}()
ep_l2 = Vector{Float64}()
ej_l2 = Vector{Float64}()
ej_hdiv = Vector{Float64}()
eφ_l2 = Vector{Float64}()

order = 2
nxs = 2:5
domain = (-0.5,0.5,-0.5,0.5,-0.5,0.5)

for n=nxs
  partition = (3,n,n)
  
  xh, trian, quad = main(
    partition=partition, order=2, domain=domain, periodic_dir=[1],
    dirichlet_tags_u=collect(1:22), dirichlet_tags_j=collect(1:22),
    f_u=f_u, f_p=f_p, f_j=f_j, f_φ=f_φ, g_u=g_u, g_j=g_j,
    ∇u_n=∇u_n, g_p=p , g_φ=φ, B=B
    )

  uh, ph, jh, φh = xh

  eu = uh - u
  ep = ph - p
  ej = jh - j
  eφ = φh - φ

  l2(v) = v*v
  h1(v) = v*v + inner(∇(v),∇(v))
  hdiv(v) = v*v + inner((∇*v),(∇*v))

  append!(eu_l2, sqrt(sum(integrate(l2(eu),trian,quad))))
  append!(eu_h1, sqrt(sum(integrate(h1(eu),trian,quad))))
  append!(ep_l2, sqrt(sum(integrate(l2(ep),trian,quad))))
  append!(ej_l2, sqrt(sum(integrate(l2(ej),trian,quad))))
  append!(ej_hdiv, sqrt(sum(integrate(hdiv(ej),trian,quad))))
  append!(eφ_l2, sqrt(sum(integrate(l2(eφ),trian,quad))))

  writevtk(trian,"results", nsubcells=10,cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh,
                                                     "eu"=>eu,"ep"=>ep,"ej"=>ej,"eφ"=>eφ])

end



slope_eu_l2 = fit([log(x) for x in nxs], [log(y) for y in eu_l2], 1)[end]
slope_eu_h1 = fit([log(x) for x in nxs], [log(y) for y in eu_h1], 1)[end]
slope_ep_l2 = fit([log(x) for x in nxs], [log(y) for y in ep_l2], 1)[end]
slope_ej_l2 = fit([log(x) for x in nxs], [log(y) for y in ej_l2], 1)[end]
slope_ej_hdiv = fit([log(x) for x in nxs], [log(y) for y in ej_hdiv], 1)[end]
slope_eφ_l2 = fit([log(x) for x in nxs], [log(y) for y in eφ_l2], 1)[end]

@test abs(slope_eu_l2 + order + 1.0) < 2e-1
@test abs(slope_eu_h1 + order) < 2e-1
@test abs(slope_ep_l2 + order) < 2e-1
@test abs(slope_ej_l2 + order + 1.0) < 2e-1
@test abs(slope_ej_hdiv + order) < 2e-1
@test abs(slope_eφ_l2 + order) < 2e-1

end #module

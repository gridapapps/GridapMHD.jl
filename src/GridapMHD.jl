module GridapMHD

using Gridap
using LineSearches: BackTracking

include("PeriodicBC.jl")
using .PeriodicBC
# using Gridap: ∇, divergence

export main
export analytical_solution

function writePVD(filename,timeSteps)
  rm(filename,force=true,recursive=true)
  mkdir(filename)
  pvdcontent  = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>\n"""
  for t in timeSteps
    pvdcontent *= """    <DataSet timestep=""" * '"'
    pvdcontent *= string(t) * '"' * """ group="" part="0" file=""" * '"'
    pvdcontent *= filename*"""/time_"""*string(t)*""".vtu"/>\n"""
  end
  pvdcontent  *= "  </Collection>\n</VTKFile>"
  f = open(filename * ".pvd", "w")
  write(f,pvdcontent)
  close(f)
end

function analytical_solution(a::Float64,       # semi-length of side walls
                             b::Float64,       # semi-length of Hartmann walls
                             t_w::Float64,     # wall thickness
                             σ_w::Float64,     # wall conductivity
                             σ::Float64,       # fluid conductivity
                             μ::Float64,       # fluid viscosity
                             grad_pz::Float64, # presure gradient
                             Ha::Float64,      # Hartmann number
                             n::Int,           # number of sumands included in Fourier series
                             x)                # evaluation point
  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  d_B = t_w*σ_w/(a*σ)

  V = 0.0
  dH_dx = 0.0; dH_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    r1_k = 0.5*( Ha + (Ha^2 + 4*α_k^2)^0.5)
    r2_k = 0.5*(-Ha + (Ha^2 + 4*α_k^2)^0.5)
    N = (Ha^2 + 4*α_k^2)^0.5

    V2 = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))+exp(-r1_k*(1+η))))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    V3 = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))+exp(-r2_k*(1+η))))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))

    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3)*(1-V2-V3)

    H2 = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))-exp(-r1_k*(1+η))))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    H3 = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))-exp(-r2_k*(1+η))))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))

    H2_dy = ((d_B * r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))) * 0.5 * (exp(-r1_k*(1-η))*(r1_k/a)-exp(-r1_k*(1+η))*(-r1_k/a)))/
         (0.5*(1+exp(-2*r1_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k)))

    H3_dy = ((d_B * r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))) * 0.5 * (exp(-r2_k*(1-η))*(r2_k/a)-exp(-r2_k*(1+η))*(-r2_k/a)))/
         (0.5*(1+exp(-2*r2_k))*d_B*N + (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k)))


    dH_dx += -2*(-1)^k*sin(α_k * ξ)/(a*l*α_k^3)*(H2-H3)
    dH_dy += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3)*(H2_dy-H3_dy)

  end
  u_z = V/μ * (-grad_pz) * a^2
  j_x = dH_dy / μ^0.5 * (-grad_pz) * a^2*σ^0.5
  j_y = -dH_dx / μ^0.5 * (-grad_pz) * a^2*σ^0.5

  u = VectorValue(0.0,0.0,u_z)
  j = VectorValue(j_x,j_y,0.0)
  return u,j
end

function main(;partition=(4,4,3),Δt=1e-4,L=1,δ=0.3,nt=4, maxit=5, exact_ic=true)


t0 = 0.0
tf = Δt*nt
domain = (-0.5,0.5,-0.5,0.5,0.0,δ)
order = 2

map(x) = VectorValue(sign(x[1])*(abs(x[1])*0.5)^0.5,   sign(x[2])*(abs(x[2])*0.5)^0.5,  x[3])
model = CartesianDiscreteModel(domain,partition,[3],map)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"dirichlet_u",append!(collect(1:21),[23,24,25,26]))
add_tag_from_tags!(labels,"dirichlet_j",append!(collect(1:20),[23,24,25,26]))

Vu = FESpace(
    reffe=:Lagrangian, order=order, valuetype=VectorValue{3,Float64},
    conformity=:H1, model=model, dirichlet_tags="dirichlet_j")

Vp = FESpace(
    reffe=:PLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)

Vj = FESpace(
    reffe=:RaviartThomas, order=order-1, valuetype=VectorValue{3,Float64},
    conformity=:Hdiv, model=model, dirichlet_tags="dirichlet_j")

Vφ = FESpace(
    reffe=:QLagrangian, order=order-1, valuetype=Float64,
    conformity=:L2, model=model, constraint=:zeromean)


trian = Triangulation(model)
degree = 2*(order+1)
quad = CellQuadrature(trian,degree)

# B? uk?
x = get_physical_coordinate(trian)
ρ = 1.0
ν = 1.0
σ = 1.0
Re = 10.0 # U = 10.0, L = 1.0/ ν = 1.0
Ha = 10.0
K = Ha / (1-0.825*Ha^(1/2)-Ha^(-1))
f(x) = VectorValue(0.0,0.0,L^3 * K / Re)
B(x) = VectorValue(0.0,Ha,0.0)
analytical_u(x) = analytical_solution(0.5,  # semi-length of side walls
                                      0.5,  # semi-length of Hartmann walls
                                      0.0,  # wall conductivity
                                      1.0,  # wall thickness
                                      1.0,  # fluid conductivity
                                      1.0,  # fluid viscosity
                                      K/Re, # presure gradient
                                      Ha,   # Hartmann number
                                      10,   # number of sumands included in Fourier series
                                      x)[1]
analytical_j(x) = analytical_solution(0.5,  # semi-length of side walls
                                      0.5,  # semi-length of Hartmann walls
                                      0.0,  # wall conductivity
                                      1.0,  # wall thickness
                                      1.0,  # fluid conductivity
                                      1.0,  # fluid viscosity
                                      K/Re, # presure gradient
                                      Ha,   # Hartmann number
                                      10,   # number of sumands included in Fourier series
                                      x)[2]

u0 = VectorValue(0.0,0.0,0.010)
gu = VectorValue(0.0,0.0,0.0)
gj = VectorValue(0.0,0.0,0.0)

U = TrialFESpace(Vu,analytical_u)
P = TrialFESpace(Vp)
j = TrialFESpace(Vj,gj)
φ = TrialFESpace(Vφ)

Y = MultiFieldFESpace([Vu, Vp, Vj, Vφ])
X = MultiFieldFESpace([U, P, j, φ])


# # xh = FEFunction(X,rand(num_free_dofs(X)))
# uh = FEFunction(U,rand(num_free_dofs(U)))
# ph = FEFunction(P,rand(num_free_dofs(P)))
# jh = FEFunction(j,rand(num_free_dofs(j)))
# φh = FEFunction(φ,rand(num_free_dofs(φ)))
#
# # uh, ph, jh, φh = xh
#
#
# # writevtk(trian,"results",cellfields=["u"=>uh, "p"=>ph, "j"=>jh, "phi"=>φh])
# writevtk(trian,"results",cellfields=["u"=>uh, "p"=>ph, "j"=>jh, "phi"=>φh])
#
# @assert false

neumanntags = [22]
btrian = BoundaryTriangulation(model,neumanntags)
degree = 2*(order-1)
bquad = CellQuadrature(btrian,degree)
nb = get_normal_vector(btrian)

if exact_ic
  un = interpolate(Vu, analytical_u)
else
  un = interpolate(Vu, u0)
end
# @law function jxB(x)
#   a = analytical_j(x)
#   b = B(x)
#   return VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])
# end
@law vprod(a,b) = VectorValue(a[2]b[3]-a[3]b[2], a[1]b[3]-a[3]b[1], a[1]b[2]-a[2]b[1])

function a(X,Y)
  u  , p  , j  , φ   = X
  v_u, v_p, v_j, v_φ = Y

  # uk*(∇(u)*v_u) + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p
  (1/Δt)*u*v_u + ν*inner(∇(u),∇(v_u)) - p*(∇*v_u) + (∇*u)*v_p - 1/ρ * vprod(j,B(x))*v_u +
  j*v_j + σ*(∇(φ)*v_j) - σ*vprod(u,B(x))*v_j - ∇(v_φ)*j
end

@law conv(u,∇u) = (∇u')*u
@law dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

c(u,v) = inner(v,conv(u,∇(u)))
dc(u,du,v) = inner(v,dconv(du,∇(du),u,∇(u)))

function l(y)
 v_u, v_p, v_j, v_φ = y
 (1/Δt)*un*v_u + v_u*f + v_p*0.0 + v_j*VectorValue(0.0,0.0,0.0) + v_φ*0.0
end

u_nbc = VectorValue(0.0,0.0,0.0)
p_nbc = 0.0
j_nbc = VectorValue(0.0,0.0,0.0)
φ_nbc = 0.0

function l_Γ(y)
 v_u, v_p, v_j, v_φ = y
 u_nbc * v_u + p_nbc * v_p + j_nbc * v_j + φ_nbc * v_φ
end

function res(X,Y)
  u   = X[1]
  v_u = Y[1]
  a(X,Y) + c(u, v_u) - l(Y)
end

function jac(X,Y,dX)
  u   = X[1]
  v_u = Y[1]
  du  = dX[1]
  a(dX,Y) + dc(u,du,v_u)
end

t_Ω = FETerm(res,jac,trian,quad)
t_Γ = FESource(l_Γ,btrian,bquad)
op  = FEOperator(X,Y,t_Ω)

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking(), iterations=maxit)
solver = FESolver(nls)

timeSteps = collect(t0:Δt:tf)

writePVD("results",timeSteps[1:end])
writevtk(trian, "results/time_"*string(t0)*".vtu",cellfields=["u"=>un])

for t in timeSteps[2:end]

  xh = solve(solver,op)
  un, pn, jn, φn = xh

  writevtk(trian,"results/time_"*string(t)*".vtu",cellfields=["u"=>un, "p"=>pn, "j"=>jn, "phi"=>φn])

  @show t

  # Update operator
  op = FEOperator(X,Y,t_Ω)
end


end
end # module

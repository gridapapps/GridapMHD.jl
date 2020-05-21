module Defaults

using Gridap

export writePVD
export shercliff_solution
export default_u_ic
export default_g
export default_∇u_n
export default_p
export default_φ
export defualt_f
export defualt_B

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

# Default force
function default_f_u(x)
  return VectorValue(0.0,0.0,0.0)
end
# Default force
function default_f_p(x)
  return 0.0
end
# Default force
function default_f_j(x)
  return VectorValue(0.0,0.0,0.0)
end
# Default force
function default_f_φ(x)
  return 0.0
end

# Default external magnetic field
function default_B(x)
  return VectorValue(0.0,1.0,0.0)
end

# Default initial conditions
function default_u_ic(x)
  return VectorValue(0.0,0.0,0.0)
end

# Default Dirichlet boundary conditions
function default_g(x)
  return VectorValue(0.0,0.0,0.0)
end

# Default Neumann boundary conditions
function default_∇u_n(x)
  return VectorValue(0.0,0.0,0.0)
end

function default_p(x)
  return 0.0
end

function default_φ(x)
  return 0.0
end


function shercliff_solution(a::Float64,       # semi-length of side walls
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




end # module

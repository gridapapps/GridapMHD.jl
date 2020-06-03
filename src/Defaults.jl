module Defaults

using Gridap

export writePVD
export shercliff_solution
export default_u_ic
export default_g
export default_p
export default_φ
export default_∇u_n
export default_f_p
export default_f_φ
export default_f_j
export defualt_f_B

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


function shercliff_u(a::Float64,       # semi-length of side walls
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

  V = 0.0; V0=0.0;
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = (Ha^2 + 4*α_k^2)^(0.5)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num1 = d_B*r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))
    num2 = (exp(-r1_k*(1-η))+exp(-r1_k*(1+η)))/2.0
    den1 = d_B*N *(1+exp(-2*r1_k))/2.0
    den2 = (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k))
    V2 = (num1 * num2)/(den1 + den2)

    num1 = d_B*r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))
    num2 = (exp(-r2_k*(1-η)) + exp(-r2_k*(1+η)))/2.0
    den1 = d_B*N*(1+exp(-2*r2_k))/2
    den2 = (1-exp(-2(r1_k+r2_k)))/(1+exp(-2*r1_k))
    V3 = (num1 * num2)/(den1 + den2)

    # V1 = 1 - N/(2*α_k^2)*((1+exp(-2*N))/(1-exp(-2*N))-exp(Ha-N)*(1+exp(-2*Ha))/(1-exp(-2*N)))

    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3) * (1-V2-V3)
    # V0+= 1/α_k^4 * V1
  end
  u_z = V/μ * (-grad_pz) * a^2
  # u_0 = -2*a^2/(l^2*η) * grad_pz * V0

  return VectorValue(0.0*u_z,0.0*u_z,u_z)
end


function shercliff_j(a::Float64,       # semi-length of side walls
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

  H_dx = 0.0; H_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = sqrt(Ha^2 + 4*α_k^2)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num1 = d_B*r2_k + (1-exp(-2*r2_k))/(1+exp(-2*r2_k))
    num2 = (exp(-r1_k*(1-η)) - exp(-r1_k*(1+η)))/2.0
    num2_dy = (exp(-r1_k*(1-η))*(r1_k/a) + exp(-r1_k*(1+η))*(r1_k/a))/2.0
    den1 = d_B*N*(1+exp(-2*r1_k))/2
    den2 = (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r2_k))
    H2 = (num1 * num2)/(den1 + den2)
    H2_dy = (num1 * num2_dy)/(den1 + den2)

    num1 = d_B*r1_k + (1-exp(-2*r1_k))/(1+exp(-2*r1_k))
    num2 = (exp(-r2_k*(1-η)) - exp(-r2_k*(1+η)))/2.0
    num2_dy = (exp(-r2_k*(1-η))*(r2_k/a) + exp(-r2_k*(1+η))*(r2_k/a))/2.0
    den1 = d_B*N*(1+exp(-2*r2_k))/2
    den2 = (1-exp(-2*(r1_k+r2_k)))/(1+exp(-2*r1_k))
    H3 = (num1 * num2)/(den1 + den2)
    H3_dy = (num1 * num2_dy)/(den1 + den2)

    H_dx += -2*(-1)^k * sin(α_k * ξ)/(a*l*α_k^2) * (H2 - H3)
    H_dy += 2*(-1)^k * cos(α_k * ξ)/(l*α_k^3) * (H2_dy - H3_dy)
  end
  j_x = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * H_dy
  j_y = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * (-H_dx)

  return VectorValue(j_x,j_y,0.0)
end




end # module

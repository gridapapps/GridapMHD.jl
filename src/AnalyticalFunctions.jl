module AnalyticalSolutions

using Gridap: VectorValue

export shercliff_u
export shercliff_j
export hunt_u
export hunt_j



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

    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3) * (1-V2-V3)
  end
  u_z = V/μ * (-grad_pz) * a^2

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


function hunt_u(a::Float64,       # semi-length of side walls
                b::Float64,       # semi-length of Hartmann walls
                μ::Float64,       # fluid viscosity
                grad_pz::Float64, # presure gradient
                Ha::Float64,      # Hartmann number
                n::Int,           # number of sumands included in Fourier series
                x)                # evaluation point
  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  V = 0.0; V0=0.0;
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = (Ha^2 + 4*α_k^2)^(0.5)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η))+exp(-r1_k*(1+η))
    den = 1+exp(-2*r1_k)
    V2 = (r2_k/N)*(num/den)

    num = exp(-r2_k*(1-η))+exp(-r2_k*(1+η))
    den = 1+exp(-2*r2_k)
    V3 = (r1_k/N)*(num/den)


    V += 2*(-1)^k*cos(α_k * ξ)/(l*α_k^3) * (1-V2-V3)
  end
  u_z = V/μ * (-grad_pz) * a^2

  return VectorValue(0.0*u_z,0.0*u_z,u_z)
end


function hunt_j(a::Float64,       # semi-length of side walls
                b::Float64,       # semi-length of Hartmann walls
                σ::Float64,       # fluid conductivity
                μ::Float64,       # fluid viscosity
                grad_pz::Float64, # presure gradient
                Ha::Float64,      # Hartmann number
                n::Int,           # number of sumands included in Fourier series
                x)                # evaluation point
  l = b/a
  ξ = x[1]/a
  η = x[2]/a

  H_dx = 0.0; H_dy = 0.0
  for k in 0:n
    α_k = (k + 0.5)*π/l
    N = sqrt(Ha^2 + 4*α_k^2)
    r1_k = 0.5*( Ha + N)
    r2_k = 0.5*(-Ha + N)

    num = exp(-r1_k*(1-η)) - exp(-r1_k*(1+η))
    num_dy = exp(-r1_k*(1-η))*(r1_k/a) + exp(-r1_k*(1+η))*(r1_k/a)
    den = 1+exp(-2*r1_k)
    H2 = (r2_k/N)*(num/den)
    H2_dy = (r2_k/N)*(num_dy/den)

    num = exp(-r2_k*(1-η)) - exp(-r2_k*(1+η))
    num_dy = exp(-r2_k*(1-η))*(r2_k/a) + exp(-r2_k*(1+η))*(r2_k/a)
    den = 1+exp(-2*r2_k)
    H3 = (r1_k/N)*(num/den)
    H3_dy = (r1_k/N)*(num_dy/den)

    H_dx += -2*(-1)^k * sin(α_k * ξ)/(a*l*α_k^2) * (H2 - H3)
    H_dy += 2*(-1)^k * cos(α_k * ξ)/(l*α_k^3) * (H2_dy - H3_dy)
  end
  j_x = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * H_dy
  j_y = a^2*σ^0.5 / μ^0.5 * (-grad_pz) * (-H_dx)

  return VectorValue(j_x,j_y,0.0)
end


end # module

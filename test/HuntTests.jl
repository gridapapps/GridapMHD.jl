module HuntTests

using GridapMHD

ρ = 1.0
ν = 1.0
σ = 1.0

U0 = 10.0
B0 = 10.0
L = 1.0

Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))

GridapMHD.hunt(nx=5,ny=5,Re=Re,Ha=Ha,debug_setup=true,title="hunt")

end # module


module MainTests

using Gridap
using GridapMHD

ρ = 1.0
ν = 1.0
σ = 1.0

U0 = 10.0
B0 = 10.0
L = 1.0
Re = U0 * L / ν
Ha = B0 * L * sqrt(σ/(ρ*ν))
N = Ha^2/Re

domain = (0,1,0,1,0,1)
cells = (10,10,10)
model = CartesianDiscreteModel(domain,cells)
Ω = Interior(model)
ud = VectorValue(0.,1.,0.)

params = Params(
  domain = Ω,
  material = Material(N=N,Re=Re),
  bcs = [
    VelocityBc(domain="boundary",value=ud),
    InsulatingBc(domain="boundary",value=ud)]
)

out = GridapMHD.main(params)

uh,ph,jh,φh = out.solution

writevtk(Ω,"test_Ω",
  nsubcells=2,
  cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh])

end # module

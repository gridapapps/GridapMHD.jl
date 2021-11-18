function hunt(;
  nx::Int=3,
  ny::Int=3,
  Re::Float64 = 10.0,
  Ha::Float64 = 10.0,
  kwargs...)

  L = 1.0
  N = Ha^2/Re
  K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1))
  ∂p∂z = -Re * K / L^3
  f_u = VectorValue(0.0,0.0, -∂p∂z) * L/Re^2 # Assumed ν=L=1 by default
  B = VectorValue(0.0,1.0,0.0)

  domain = (-1.0,1.0,-1.0,1.0,0.0,0.1)
  map(x) = (2/sqrt(2))*VectorValue(
    sign(x[1])*(abs(x[1])*0.5)^0.5,
    sign(x[2])*(abs(x[2])*0.5)^0.5,
    x[3]*sqrt(2)/2)

  partition=(nx,ny,3)
  model = CartesianDiscreteModel(
    domain,partition;isperiodic=(false,false,true),map=map)
  Ω = Interior(model)

  labels = get_face_labeling(model)
  tags_u = append!(collect(1:20),[23,24,25,26])
  tags_j = append!(collect(1:20),[25,26])
  tags_φ = [23,24]
  add_tag_from_tags!(labels,"noslip",tags_u)
  add_tag_from_tags!(labels,"insulating",tags_j)
  add_tag_from_tags!(labels,"conducting",tags_φ)

  actions = [
    ConductingFluid(domain=Ω,α=1,β=(1/Re)),
    VelocityBc(domain="noslip"),
    InsulatingBc(domain="insulating"),
    ConductingBc(domain="conducting"),
    MagneticField(domain=Ω,B=B,γ=N),
    ]

  out = main(model,actions;kwargs...)

  out
end

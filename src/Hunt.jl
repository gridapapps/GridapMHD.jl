function hunt(;
  nx::Int=3,
  ny::Int=3,
  Re::Float64 = 10.0,
  Ha::Float64 = 10.0,
  title::String="test",
  vtk::Bool=false)

  L = 1.0
  N = Ha^2/Re
  K = Ha / (1-0.95598*Ha^(-1/2)-Ha^(-1))
  ∂p∂z = -Re * K / L^3
  f_u = VectorValue(0.0,0.0, -∂p∂z) * L/Re^2 # Assumed ν=L=1 by default
  B0 = VectorValue(0.0,1.0,0.0)

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
  dirichlet_tags_u = append!(collect(1:20),[23,24,25,26])
  dirichlet_tags_j = append!(collect(1:20),[25,26])
  add_tag_from_tags!(labels,"dirichlet_u",dirichlet_tags_u)
  add_tag_from_tags!(labels,"dirichlet_j",dirichlet_tags_j)

  bc1 = VelocityBc(domain="dirichlet_u")
  bc2 = InsulatingBc(domain="dirichlet_j")

  params = Params(
    domain = Ω,
    material=Material(N=N,Re=Re),
    bcs = [bc1,bc2],
    B0 = B0,
  )

  out = main(params)

  if vtk
    writevtk(Ω,"$(title)_Ω",
      nsubcells=2,
      cellfields=["uh"=>uh,"ph"=>ph,"jh"=>jh,"φh"=>φh])
  end

  out
end

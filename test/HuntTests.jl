module HuntTests

using GridapMHD
using Gridap

info = GridapMHD.hunt(
  nx=10,
  ny=10,
  L=1.0,
  B=VectorValue(0.,50.,0.),
  debug=false,
  vtk=true,
  title="hunt")
display(info)

end # module


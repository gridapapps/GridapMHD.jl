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
  title="hunt",
  solver="petsc",
  petsc_options="-snes_monitor -ksp_monitor -ksp_view -pc_type lu"
 )
display(info)

end # module


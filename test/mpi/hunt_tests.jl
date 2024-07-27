module HuntTestsMPI

using GridapPETSc
using SparseMatricesCSR

using GridapMHD: hunt

function main(parts)
  # Default monolithic solver w petsc
  hunt(
    nc=(4,4),
    np=parts,
    backend=:mpi,
    L=1.0,
    B=(0.,50.,0.),
    debug=false,
    vtk=true,
    title="hunt",
    solver=:petsc,
  )
end

end
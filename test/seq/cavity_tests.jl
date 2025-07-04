module CavityTestsSequential

using GridapPETSc
using SparseMatricesCSR, SparseArrays
using GridapMHD: cavity

# Serial, LUSolver
cavity(
  fluid_disc = :Qk_dPkm1,
  current_disc = :RT,
)

cavity(
  fluid_disc = :Qk_dPkm1,
  current_disc = :H1,
)

cavity(
  fluid_disc = :RT,
  current_disc = :RT,
)

end # module

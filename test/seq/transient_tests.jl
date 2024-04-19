module TransientTestsSequential

using GridapMHD: transient
using GridapPETSc, SparseMatricesCSR

transient(;man_solution=:exact,Δt=0.1,tf=0.1)

end # module

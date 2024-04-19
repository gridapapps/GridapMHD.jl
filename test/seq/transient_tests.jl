module TransientTestsSequential

using GridapMHD: transient
using GridapPETSc, SparseMatricesCSR

transient(;man_solution=:exact,Δt=0.1,tf=0.1,max_error=1e-8)

end # module

module TransientTestsSequential

using GridapMHD: transient
using GridapPETSc, SparseMatricesCSR

transient(;man_solution=:stationary_fespace,Δt=0.1,tf=0.1,max_error=1e-8)
transient(;man_solution=:lineartime_fespace,Δt=0.1,tf=0.1,max_error=1e-8)

end # module

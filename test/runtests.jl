
using Test

@info "Sequential tests"
include("seq/runtests.jl")

@info "MPI tests"
#MPI.mpiexec() do cmd
#  cmd = `$(cmd) -np 4 julia -e 'include("test/mpi/runtests.jl")'`
#end

module ExpansionTestsSequential

using GridapMHD: expansion

expansion()
expansion(np=2,backend=:sequential)
expansion(np=1,backend=:mpi)

end # module

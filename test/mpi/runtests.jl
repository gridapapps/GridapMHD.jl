module GridapMHDMPITests

using Test
using MPI

testdir = @__DIR__
istest(f) = endswith(f, ".jl") && !(f=="runtests.jl")
testfiles = sort(filter(istest, readdir(testdir)))

MPI.mpiexec() do cmd
  for file in testfiles
    path = joinpath(testdir,file)
    _cmd = `$(cmd) -np 4 julia --project=. $path`
    run(_cmd)
  end
end

end
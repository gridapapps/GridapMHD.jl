module GridapMHDMPITests

using Test
using MPI

tests_to_run = ["cavity_tests.jl","hunt_tests.jl","expansion_tests.jl"]

testdir = @__DIR__
istest(f) = endswith(f, ".jl") && !(f=="runtests.jl") && (f in tests_to_run)
testfiles = sort(filter(istest, readdir(testdir)))

MPI.mpiexec() do cmd
  for file in testfiles
    path = joinpath(testdir,file)
    _cmd = `$(cmd) -np 4 --allow-run-as-root --oversubscribe $(Base.julia_cmd()) --project=. $path`
    @show cmd
    run(_cmd)
    @test true
  end
end

end
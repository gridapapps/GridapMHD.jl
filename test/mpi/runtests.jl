module GridapMHDMPITests

using Test
using MPI


testdir = @__DIR__
# istest(f) = endswith(f, ".jl") && !(f == "runtests.jl")
# testfiles = sort(filter(istest, readdir(testdir)))

function run_test(n, file)
  MPI.mpiexec() do cmd
      path = joinpath(testdir, file)
      _cmd = `$(cmd) -np $n --allow-run-as-root --oversubscribe $(Base.julia_cmd()) --project=. $path`
      @show _cmd
      run(_cmd)
      @test true
  end
end

run_test(1,"runtests_body.jl")
# run_test(4,"runtests_body.jl")

end
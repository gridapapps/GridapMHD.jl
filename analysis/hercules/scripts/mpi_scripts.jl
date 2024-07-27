using MPI

projectdir = joinpath(@__DIR__,"..","..","..") |> abspath
sysimage = joinpath(projectdir,"GridapMHD.so")

function num_procs(;np,backend,kwargs...)
  @assert backend == :mpi
  prod(np)
end

function taskcmd(;app="channel",precompile=false,setup="",solver=":julia",kwargs...)
  kw = str_kwargs(;kwargs...)
  prec = precomp(;precompile,kwargs...)
  str = "
  using MPI
  MPI.Init()
  if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println(\"Num Procs:\",MPI.Comm_size(MPI.COMM_WORLD))
  end
  using GridapMHD: $app
  $setup
  solver = $solver
  $(prec)
  @time $app(;$kw,solver=solver)
  "
  str
end

function precomp(;precompile,app="channel",np=(1,1,1),nc=np,solver=:julia,kwargs...)
  str = ""
  if precompile
    nc=np.*2
    str *= "@time $app(;np=$np,nc=$nc,solver=:$solver,backend=:mpi)"
  end
  str
end

function str_kwargs(;tabsize=4,kwargs...)
  NamedTuple(kwargs)
  tab = ' ' ^ tabsize
  str = string(NamedTuple(kwargs))
  str = str[2:end-1]
  str
end

function run_driver(;kwargs...)
  args = `-O3 --check-bounds=no -J$sysimage --project=.`
  args = `-O3 --check-bounds=no --project=.`
  args = `--project=$projectdir`
  MPI.mpiexec() do cmd
    n = num_procs(;kwargs...)
    task = taskcmd(;kwargs...)
    _cmd = `$(cmd) -np $n --oversubscribe $(Base.julia_cmd()) $args -e $task`
    @show _cmd
    run(_cmd)
  end
end

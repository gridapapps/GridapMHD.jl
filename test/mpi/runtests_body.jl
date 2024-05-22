using MPI

if ! MPI.Initialized()
  MPI.Init()
end
if MPI.Comm_size(MPI.COMM_WORLD) == 4
  parts=(2,2)
elseif MPI.Comm_size(MPI.COMM_WORLD) == 1
  parts=(1,1)
else
  MPI.Abort(MPI.COMM_WORLD,0)
end

# include("hunt_tests.jl")
# HuntTestsMPI.main(parts)

include("hunt_li2019_tests.jl")
HuntLi2019TestsMPI.main(parts)

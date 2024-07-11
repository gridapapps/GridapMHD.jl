module ChannelTestsSequential

using GridapMHD: channel
using GridapPETSc, SparseMatricesCSR


channel(nc=(4,4,4),vtk=true)

channel(backend=:mpi,np=(1,1,1),solver=:petsc,
  sizes=(4,2,2),nc=(4,4,4),vtk=true)

channel(backend=:mpi,np=(1,1,1),solver=:petsc,
  sizes=(8,2,2),nc=(4,8,8),B0=20,vtk=true,inlet=:shercliff)

channel(sizes=(4,1,1),nc=(16,4,4),vtk=true)

channel(sizes=(4,1,1),nc=(16,8,8,),debug=false,vtk=true,
  bl_orders=(1,2,2),nonuniform_B=true,B0=20,inlet=:shercliff)


#TODO: be careful with diverging computations


# channel(
#   backend=:mpi,np=(1,1,1),solver=:petsc,
#   sizes=(8,2,2),nc=(10,10,10),vtk=true,inlet=:shercliff,bl_orders=(1,2,2),
#   B0=100/2,u0=10/2)

# channel(
#   title="osc_1",
#   backend=:mpi,np=(1,1,1),solver=:petsc,
#   sizes=(8,2,2),nc=(15,15,15),vtk=true,inlet=:shercliff,bl_orders=(1,2,2),
#   B0=50/2,u0=10/2,niter=1)


#TODO: rename channel to channel

# channel(backend=:mpi,np=(1,1,1),sizes=(4,1,1),nc=(8,8,8),ν=1,vtk=true,
# B0=10/4,bl_orders=(1,2,2),μ=1e-2)

# nc = (16,10,10) # Limit of HERCULES
# nc = (8,8,8)
# L = 4
# γ = 0.45*30/L
# channel(backend=:mpi,np=(1,1,1),sizes=(L,1,1),nc=nc,ν=1,vtk=true,
# bl_orders=(1,2,2),u0=10/L,B0=100/L,nonuniform_B=true,γB=γ,μ=1e-2)


# channel(backend=:mpi,np=(1,1,1),sizes=(4,1,1),nc=(16,4,4),ν=1e-1,vtk=true,
# nonuniform_B=true)

# ## Stabilization
# L = 4
# n = 6
# # Re = 100, Ha = 1000, μ = 0, niter= 1 (not stable)
# channel(
#   backend=:mpi,
#   np=(1,1,1),
#   sizes=(L,1,1),
#   nc=(L*n,n,n),
#   ν=1,
#   B0=1000/L,
#   u0=50/L,
#   μ=0,
#   niter=1,
#   vtk=true)
# # Re = 100, Ha = 1000, μ = 2e-2, niter= 1 (stable)
# channel(
#   backend=:mpi,
#   np=(1,1,1),
#   sizes=(L,1,1),
#   nc=(L*n,n,n),
#   ν=1,
#   B0=1000/L,
#   u0=50/L,
#   μ=100,
#   niter=1,
#   vtk=true)



# L = 4
# n = 8
# channel(
#   backend=:mpi,
#   np=(1,1,1),
#   sizes=(L,1,1),
#   nc=(L*n,n,n),
#   ν=1e-2,
#   B0=100/L,
#   u0=1/L,
#   μ=2e-2,
#   nonuniform_B=true,
#   vtk=true)

end

module PipeTestsSequential

using GridapMHD: pipe
using GridapPETSc, SparseMatricesCSR


pipe(vtk=true)


pipe(backend=:mpi,np=(1,1,1),solver=:petsc,
  sizes=(4,1,1),nc=(4,4,4),ν=1e-1,vtk=true)

pipe(sizes=(4,1,1),nc=(16,4,4),ν=1e-1,vtk=true)

pipe(sizes=(4,1,1),nc=(16,4,4),ν=1e-3,vtk=true)
pipe(sizes=(4,1,1),nc=(16,4,4),ν=1e-2,μ=2e-2,vtk=true)
pipe(sizes=(4,1,1),nc=(16,4,4),ν=5e-3,μ=1e-1,vtk=true)

pipe(sizes=(4,1,1),nc=(16,8,8,),debug=true,vtk=true,
  bl_orders=(1,2,2),nonuniform_B=true)


# pipe(backend=:mpi,np=(1,1,1),sizes=(4,1,1),nc=(8,8,8),ν=1,vtk=true,
# B0=10/4,bl_orders=(1,2,2),μ=1e-2)

# nc = (16,10,10) # Limit of HERCULES
# nc = (8,8,8)
# L = 4
# γ = 0.45*30/L
# pipe(backend=:mpi,np=(1,1,1),sizes=(L,1,1),nc=nc,ν=1,vtk=true,
# bl_orders=(1,2,2),u0=10/L,B0=100/L,nonuniform_B=true,γB=γ,μ=1e-2)


# pipe(backend=:mpi,np=(1,1,1),sizes=(4,1,1),nc=(16,4,4),ν=1e-1,vtk=true,
# nonuniform_B=true)

# ## Stabilization
# L = 4
# n = 6
# # Re = 100, Ha = 1000, μ = 0, niter= 1 (not stable)
# pipe(
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
# pipe(
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
# pipe(
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

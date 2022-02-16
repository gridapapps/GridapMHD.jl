module load openmpi/4.1.2
module load intel-mkl/2022.0.2
module load julia/1.7.1
export PETSCROOT=/home/552/fv3851/Apps/petsc/3.15.5
export JULIA_PETSC_LIBRARY=/home/552/fv3851/Apps/petsc/3.15.5/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/apps/openmpi/4.1.2
export GRIDAPMHD=/home/552/fv3851/Code/GridapMHD.jl


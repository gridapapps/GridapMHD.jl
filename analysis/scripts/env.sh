module load openmpi/4.1.2
module load intel-mkl/2022.0.2
module load julia/1.7.1
export PETSCROOT=/home/552/fv3851/Apps/petsc/3.15.5
export JULIA_PETSC_LIBRARY=/home/552/fv3851/Apps/petsc/3.15.5/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/apps/openmpi/4.1.2
export GRIDAPMHD=/scratch/bt62/fv3851/GridapMHD.jl
export UCX_WARN_UNUSED_ENV_VARS=n
export HCOLL_ML_DISABLE_SCATTERV=1
export HCOLL_ML_DISABLE_BCAST=1

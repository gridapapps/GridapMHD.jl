module load openmpi_gcc/4.0.4

export PETSCROOT=/home/froca/software/petsc/3.15.5
export JULIA_PETSC_LIBRARY=/home/froca/software/petsc/3.15.5/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/opt/openmpi_gcc/4.0.4
export UCX_WARN_UNUSED_ENV_VARS=n
export HCOLL_ML_DISABLE_SCATTERV=1
export HCOLL_ML_DISABLE_BCAST=1

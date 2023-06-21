module load openmpi/4.1.2
module load intel-mkl/2022.0.2
module load julia
#export PATH=${HOME}/bin:$PATH # add (link to) julia executable here
export PETSCROOT=${HOME}petsc-3.18.4
export JULIA_PETSC_LIBRARY=${HOME}petsc-3.18.4/gcc-8.5.0_opmpi-4.1.2_mkl-2022.0.2_opt/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/apps/openmpi/4.1.2
export GRIDAPMHD=${HOME}/GridapMHD.jl
export UCX_WARN_UNUSED_ENV_VARS=n
export HCOLL_ML_DISABLE_SCATTERV=1
export HCOLL_ML_DISABLE_BCAST=1

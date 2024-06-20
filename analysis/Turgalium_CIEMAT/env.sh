module load openmpi_gcc/4.0.4
module load gcc/11.2.0


export PETSCROOT=/home/froca/software/petsc/3.19.3
export JULIA_PETSC_LIBRARY=/home/froca/software/petsc/3.19.3/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/apps/openmpi_mlnx/4.0.4
export GMSHROOT=/home/froca/software/gmsh/4.10.5
export GRIDAPMHD=/ws/blankets/GridapMHD.jl

export OMPI_MCA_btl='^openib'
export OMPI_MCA_opal_warn_on_missing_libcuda=0

#export UCX_WARN_UNUSED_ENV_VARS=n
#export HCOLL_ML_DISABLE_SCATTERV=1
#export HCOLL_ML_DISABLE_BCAST=1

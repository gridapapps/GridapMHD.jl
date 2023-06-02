module load openmpi_gcc/4.0.4
module load gcc/11.2.0


export PETSCROOT=/home/froca/software/petsc/3.18.0_new
export JULIA_PETSC_LIBRARY=/home/froca/software/petsc/3.18.0_new/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/opt/openmpi_gcc/4.0.4
export GMSHROOT=/home/froca/software/gmsh/4.10.5
export GRIDAPMHD=/ws/blankets/GridapMHD.jl

export OMPI_MCA_btl='^openib'



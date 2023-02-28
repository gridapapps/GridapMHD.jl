module load openmpi_gcc/4.0.4

export PETSCROOT=/home/froca/software/petsc/3.18.0
export JULIA_PETSC_LIBRARY=/home/froca/software/petsc/3.18.0/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/opt/openmpi_gcc/4.0.4/
export GMSHROOT=/home/froca/software/gmsh/4.10.5
export GRIDAPMHD=/home/froca/blankets/GridapMHD.jl

export OMPI_MCA_btl='^openib'



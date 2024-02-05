module load gcc-13.1.0-gcc-8.5.0-ueufrru
module load openmpi-4.1.5-gcc-13.1.0-nzcz542 

export PETSCROOT=/mnt/lustre/home/u6678/software/petsc/3.20.2
export JULIA_PETSC_LIBRARY=/mnt/lustre/home/u6678/software/petsc/3.20.2/lib/libpetsc.so
export JULIA_MPI_BINARY=system
export JULIA_MPI_PATH=/mnt/lustre/home/spack/spack/opt/spack/linux-rocky8-icelake/gcc-13.1.0/openmpi-4.1.5-nzcz542rea5aajaza2x6ps3ym4krutku/
export GMSHROOT=/mnt/lustre/home/u6678/software/gmsh/4.11.1
export GRIDAPMHD=/mnt/lustre/home/u6678/work/GridapMHD.jl

#export PMIX_MCA_psec=^munge

#export OMPI_MCA_btl='^openib'
#export OMPI_MCA_opal_warn_on_missing_libcuda=0

#export UCX_WARN_UNUSED_ENV_VARS=n
#export HCOLL_ML_DISABLE_SCATTERV=1
#export HCOLL_ML_DISABLE_BCAST=1

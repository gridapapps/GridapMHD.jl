#!/bin/sh

#SBATCH -N 1 
#SBATCH --ntasks-per-node=2
#SBATCH -t 1:00:00
#SBATCH --partition=compile

#SBATCH -o output.txt 
#SBATCH -e error.txt 
#SBATCH --job-name=compilation
#SBATCH --chdir=/home/froca/blankets/GridapMHD.jl/analysis

SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs"
echo "================================================================"
echo ""

export MALLOC_CHECK_=1

SECONDS=0

source ./env.sh

julia --project=.. -e 'using Pkg; Pkg.build(verbose=true)'

mpiexec -n 1 julia --project=.. -O3 --check-bounds=no --color=yes ../compile/compile.jl

duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID
STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""

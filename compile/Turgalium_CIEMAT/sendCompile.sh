#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 20:00:00
#SBATCH --partition=volta

#SBATCH -o outputCompile36c
#SBATCH -e errorCompile36c
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=compileGridapMHD36c
#SBATCH --mem=0


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs in ${SLURM_JOB_NUM_NODES} nodes"
echo "================================================================"
echo ""


SECONDS=0


source /ws/blankets/GridapMHD.jl/analysis/Turgalium_CIEMAT/env.sh

#julia --project=$GRIDAPMHD -e 'using Pkg; Pkg.precompile()'
julia --project=$GRIDAPMHD -O3 --check-bounds=no --color=yes $GRIDAPMHD/compile/Turgalium_CIEMAT/compile36c.jl

duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""


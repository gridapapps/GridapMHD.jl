#!/bin/bash
#PBS -q normal 
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=190gb
#PBS -l jobfs=1gb
#PBS -N compile
#PBS -l wd

source env.sh

julia --project=$GRIDAPMHD -O3 --check-bounds=no --color=yes $GRIDAPMHD/compile/compile.jl


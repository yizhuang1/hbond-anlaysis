#!/bin/bash
#SBATCH --job-name="hb_tele4"
#SBATCH --output="namd.%j.%N.out"
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=24:00:00


# job_________________________

#module load python/2.7


# run job
python weighthb.py
python weighthb2.py

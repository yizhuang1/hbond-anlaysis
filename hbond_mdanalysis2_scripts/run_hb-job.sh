#!/bin/bash
#SBATCH --job-name="ala30_hb"
#SBATCH --output="namd.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --no-requeue
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH -t 48:00:00

# job_________________________


module load python/3.8

# run job
python2 run_hbpkl_post.py                         

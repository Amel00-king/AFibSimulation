#!/bin/bash
#SBATCH -p compute          # Submit to 'compute' Partitiion or queue\
#SBATCH -J Electrode_sim          # Name the job as 'MPItest'\
#SBATCH -o Electrode_sim-%j.out   # Write the standard output to file named 'jMPItest-<job_number>.out'\
#SBATCH -e Electrode_sim-%j.err   # Write the standard error to file named 'jMPItest-<job_number>.err'\
#SBATCH -t 4-15:00:00        # Run for a maximum time of 4 days, 15 hours, 00 mins, 00 secs\
#SBATCH --nodes=1            # Request N nodes\
#SBATCH --ntasks-per-node=8 # Request n cores or task per node\

module load Miniconda   # Load the compiler and  the MPI library \
module list                 # will list modules loaded; we'll just use this to check that the modules we selected are indeed loaded\
pwd                         # prints current working directory\
date                        # prints the date and time\
python3 ./numba_tom_backwards.py    # run the MPI job\

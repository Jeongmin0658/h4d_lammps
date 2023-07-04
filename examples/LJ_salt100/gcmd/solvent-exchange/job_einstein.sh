#!/bin/bash
# temps CPU a ajuster au calcul
#PBS -l walltime=48:00:00
# memoire a ajuster au calcul
##PBS -l vmem=5gb
# nb procs a ajuster au calcul, voir OMP_NUM_THREADS plus bas
#PBS -l nodes=1:ppn=27
# job name - a changer, nom du job (qstat)
#PBS -N lj_test
#Request that regular output and terminal output go to the same file
#PBS -j oe
#  Request that your login shell variables be available to the job
#PBS -V

# einstein
cd $PBS_O_WORKDIR
# ------------ tmp running --------------------
module purge
module load lammps/8Feb2023

export INTEL="-pk intel 0 omp 1 mode mixed -sf intel"
export KMP_BLOCKTIME=0  # for intel package
#export ACCEL_CMD="-pk omp 2 -sf omp"

# Pour savoir sur quel noeud on est
echo "Working dirctory	: $PBS_O_WORKDIR"

cpu="27"
# lammps + python
run_file="run_gcmc_sample.py" #"run_hneMDMC_lammps.py"
#mypython="/SSD/jmkim/Program/anaconda3/bin/python3.10"
mypython="/SSD/jmkim/Program/anaconda3/envs/py3.7/bin/python"

date
#####################################
## run 
pwd
mpirun ${mypython} ${run_file} # $INTEL
#####################################
date

#mv ~/*.o* .

#exit


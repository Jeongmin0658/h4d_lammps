#!/bin/bash
# temps CPU a ajuster au calcul
#PBS -l walltime=10:00:00
# memoire a ajuster au calcul
##PBS -l vmem=5gb
# nb procs a ajuster au calcul, voir OMP_NUM_THREADS plus bas
#PBS -l nodes=1:ppn=8
# job name - a changer, nom du job (qstat)
#PBS -N LJ_mu_salt_deletion
#Request that regular output and terminal output go to the same file
#PBS -j oe
#  Request that your login shell variables be available to the job
#PBS -V

cd $PBS_O_WORKDIR

# Pour savoir sur quel noeud on est
echo "Working dirctory	: $PBS_O_WORKDIR"

cpu="8"
# lammps + python
run_file="run_gcmc_sample.py" #"run_hneMDMC_lammps.py"
mypython="/SSD/jmkim/Program/anaconda3/envs/py3.7/bin/python"

date
#####################################
## run 
pwd
mpirun ${mypython} ${run_file}
#####################################
date

exit


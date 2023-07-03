#############################
# Basic modules
#############################
#from __future__ import print_function
from ctypes import *        
import sys,random,math
import numpy as np
import os

# import dictionary from this directory
# dictionary for inputs and data
from dictionaries import LMP_dics
mydics=LMP_dics()

#sys.path.append('/data4/jeongmin/Simulations/separate_python_sources_ver2')

#############################
# MPI setting
#############################
me = 0
from mpi4py import MPI
me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

#############################
# Load LAMMPS
#############################
from lammps import lammps

#############################
# Seed random number generator
#############################
rseed=210900
random.seed(rseed)

#############################
# Import modules
#############################
## dictionary for inputs and data
#from dictionaries import LMP_dics
#mydics=LMP_dics()
# misc. subroutines
from define_lmp_misc import LMP_misc
mymisc=LMP_misc()
# define force field
from define_pair import LMP_pair
mypair=LMP_pair()
# MD propagation driver
from define_propagate import LMP_propagate
mymd=LMP_propagate()
# MC trial moves
from define_trial_moves import LMP_trial_moves
mymove=LMP_trial_moves()
# GCMC particle type
from define_flying import LMP_flying
myflying=LMP_flying()
# MC driver
from define_MC import LMP_MC
mymonte=LMP_MC()
# main driver
from control_ver7 impot LMP_control
#mycontrol=LMP_control(rseed=rseed)

if __name__=="__main__":
        # my driver for hneMDMC
        mycontrol=LMP_control(rseed=rseed)
        # run hybrid neMD/MC
        mycontrol.run_hneMDMC(
            read_input=True,
            n_preMD=10000,
            n_MC=10000,
            n_eqMD=10000,
            n_neqMD=12000, #60000,
            nwaters=80, #100,
            nsalts=0,
            region_to_insert="confined",
            testprint=True,
            filemc="mc_statistics.out")


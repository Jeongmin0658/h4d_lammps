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

#############################################################
### CHANGE HERE ###
# where your source folder is
source_where="WHERE_YOUR_SOURCES_ARE"
sys.path.append(source_where)
#############################################################

###############
# main driver
###############
from control import LMP_control
# my driver for hneMDMC
mycontrol=LMP_control(cut_sel=100, width_sel=1) #rseed=rseed)

if __name__=="__main__":
        mycontrol.run_hneMDMC(
            save_freq=1,                # saving H4D outputs
            n_preMD=0,                  # number of steps in MD before GCMD
            n_MC_s=0,                   # no use this
            n_MC=10,                    # total number of MC attemps
            n_MC_each=1,                # number of particle-exchange attempts in a single MC step (sampling purpose)
            n_eqMD=200,                 # number of steps in equilibrium MD
            n_neqMD=200,                # number of stpes in NEMD
            region_to_insert="NULL",    # Do NOT change this
            testprint=False,            # Do NOT change this
            filemc="mc_statistics.out") # output file

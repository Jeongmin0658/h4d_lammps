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

#############################
# where the source folder is
#############################
# source folder
source_where="WHERE_YOUR_SOURCES_ARE"
sys.path.append(source_where)

###############
# main driver
###############
from control import LMP_control
# my driver for hneMDMC
mycontrol=LMP_control(cut_sel=50,width_sel=0.1, selection=False)

if __name__=="__main__":
        mycontrol.run_hneMDMC(
            read_input=True,
            save_freq=2,
            n_preMD=0,
            n_MC_s=0,
            n_MC=10,
            n_MC_each=1,
            n_eqMD=1000,   # to calculate eq PES
            n_neqMD=5000, #60000,
            region_to_insert="NULL",
            filemc="mc_statistics.out")

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
# anubis
source_where="/SSD/jmkim/Simulation/H4D_method/Github/Github_h4d_Ongoing/ver0.0"
sys.path.append(source_where)

###############
# main driver
###############
from control import LMP_control
# my driver for hneMDMC
mycontrol=LMP_control(cut_sel=50,width_sel=0.1,selection=False, rseed=1009)

if __name__=="__main__":
        mycontrol.run_hneMDMC(
            read_input=True,
            save_freq=100,
            n_preMD=0,
            n_MC_s=0,
            n_MC=5,
            n_MC_each=5,
            n_eqMD=2000,   # to calculate eq PES
            n_neqMD=10000, #10000, #60000,
            region_to_insert="NULL",
            testprint=True,
            filemc="mc_statistics.out")

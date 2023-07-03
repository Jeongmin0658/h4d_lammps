#############################
# Basic modules
#############################
#from __future__ import print_function
from ctypes import *        
import sys,random,math
import numpy as np
import os

###############################
### MPI setting
###############################
#me = 0
#from mpi4py import MPI
#me = MPI.COMM_WORLD.Get_rank()
#nprocs = MPI.COMM_WORLD.Get_size()
#
#############################
# Load LAMMPS
#############################
from lammps import lammps, MPIAbortException

#############################
# Seed random number generator
#############################
#rseed=210900
#random.seed(rseed)

#############################
# Import modules
#############################
# dictionary for inputs and data
from dictionaries import LMP_dics
mydics=LMP_dics()

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
# Bias salt
from define_bias_salt import LMP_define_bias_salt
mysalt=LMP_define_bias_salt()
# GCMC particle type
from define_flying import LMP_flying
myflying=LMP_flying()
# MC driver
from define_MC import LMP_MC
mymonte=LMP_MC()

#############################
# Control driver for hybrid neMD/MC
#############################
class LMP_control():
    def __init__(self, cut_sel=12, width_sel=0.1, selection=False, rseed=101010, lammps_name=""):
        self.intro="Hybrid NEMD/MC method\nVersion: ver0.0\nDate: 29/July/2022\nAuthor: Jeongmin Kim @ PHENIX, CNRS"
        # prep selection
        mymonte.cut_sel=cut_sel
        mymonte.width_sel=width_sel
        mymonte.selection=selection
        # random numer
        random.seed(rseed)

        # Initialize LAMMPS
        lmp=lammps(name=lammps_name) #,cmdargs=["-log","none","-screen","none"])
        lmp.version()

        # make it accessible
        self.lmp=lmp
        # computes
        mymisc.calc_thermo(mylmp=lmp)

        # Save input dictionaries
        self.system_setup(mylmp=lmp)

    def system_setup(self, mylmp=None):
        # update system setup
        self.dic_mc=mydics.dic_mc
        self.dic_flying=mydics.dic_flying
        self.dic_inputs=mydics.dic_inputs
        self.dic_peratom=mydics.dic_peratom
        # list of dictionaries to keep track of their updated values
        self.dic_list={
            "mc": self.dic_mc,
            "flying": self.dic_flying,
            "inputs": self.dic_inputs,
            "peratom": self.dic_peratom
        }   #[self.dic_mc,self.dic_flying,self.dic_inputs,self.dic_peratom]
        if not self.dic_mc['lj_flag']:
            self.thermo=mymisc.kT_cal/298.*self.dic_inputs['temp']
            #self.vol_factor=1.
        else:   # lj
            self.thermo=self.dic_inputs['temp'] # kb = 1
        # in case of confinement with a charged wall
        self.vol_factor=self.dic_inputs['confinement_factor']
        # thermal energy kbT
        self.dic_inputs['thermal']=self.thermo  # kbT

        # update random seed
        random.seed(self.dic_inputs["random_seed"])
        mymisc.print_out(mylmp=mylmp,m="system setup is updated")

    def read_input(self,mylmp=None,output="mc_test"):
        # read input file
        file_where=self.dic_inputs['input_file']['where']
        file_to_read=self.dic_inputs['input_file']['to_read']
        filelmp=file_where+"/"+file_to_read
        mylmp.file(filelmp)
        self.gcmc_input_update(mylmp=mylmp,output=output)

    def gcmc_input_update(self,mylmp=None,output="mc_test"):
        print ("=========================================",file=output)
        # update volume
        volume=mymisc.get_volume(mylmp=mylmp)
        self.dic_mc['volume']=volume*self.vol_factor # in cubic anstrom
        print ("Volume in cubic A (factor=%f):"%self.vol_factor,mymisc.volume,self.dic_mc['volume'],file=output)
        # activity update
        list_mol=['water','salt','lj']
        for item in list_mol:
            act=self.dic_flying[item]['activity']
            self.dic_mc['activity'][item]=act #math.exp(act) #/self.thermo)
            print ("Activity (chemical potential of %s):"%item,item,self.dic_mc['activity'][item],file=output)
        print ("=========================================",file=output)

    def check_if_instant_MC(self, n_neqMD=0, flag_instant=False, output=None):
        if self.dic_mc['select_salt']: 
            # either instant MC vs NEMD/MC
            if self.dic_mc['select_insert']:
                if flag_instant:
                    n_neqMD_sel=0
                else:
                    n_neqMD_sel=n_neqMD
            else:
                # need to fix this part according to a bias function in use
                p_instant=mysalt.prob_bias_select_instant(np.power(self.dic_list['mc']['volume'],1./3.),self.dic_mc['bias']['alpha'],self.dic_mc['bias']['ex_vol'], x_center=self.dic_mc['bias']["center"], bias_f=self.dic_mc['bias']["bias_f"])
                print ("excluded volume for instataneous deletion: %f"%p_instant,file=output)
                if random.random() < p_instant:
                    n_neqMD_sel=0
                    flag_instant=True
                else:
                    n_neqMD_sel=n_neqMD
                    flag_instant=False
            mydics.update_noneq_altitude_protocol(nstep=n_neqMD_sel,dic_inputs=self.dic_inputs,output=output)
        else:   # no instant MC for water or lj particle exchange
            n_neqMD_sel=n_neqMD
        return n_neqMD_sel,flag_instant

    def run_hneMDMC(self,
        read_input=True,            # input LAMMPS file to read
        zero_salt_del=False,        # only for the case calculating the salt chemical potential at infinite dilution
        save_freq=10,               # frequency to save restart files
        n_MC_s=0,                   # starting index for the exchange MC
        n_MC=10,                    # Number of MC for the exchange
        n_MC_each=1,                # Number of the exchange for a single configuration; this could be more than one when calculating the chemical potential 
        n_preMD=10000,              # Number of steps for MD before actual GCMD
        n_eqMD=1000,                # Number of steps in Equilibrium phase
        n_neqMD=12000,              # Number of steps in NEq phase for the exchange
        testprint=True,             # print option
        region_to_insert="NULL",    # Defining the region for a trial insertion move
        filemc="mc_statistics.out"  # file name for the information of GCMD
        ):
        # my lammps
        mylmp=self.lmp
        # open file for MC statistics
        mcout=open(filemc,'w')
        print (self.intro, file=mcout)
        # read input
        if read_input:
            self.read_input(mylmp=mylmp,output=mcout)
        else:
            self.gcmc_input_update(mylmp=mylmp,output=mcout)
        # additional MD run
        mymd.run_MD(mylmp=mylmp,ff=True,md="eqMD-pre",nsteps=n_preMD,dic_list=self.dic_list,output=mcout,mpiexcept=MPIAbortException)
        # update non-eq insertion/deletion velocity
        mydics.update_noneq_altitude_protocol(nstep=n_neqMD,dic_inputs=self.dic_inputs,output=mcout)
        # set list of existing water/salt/lj
        myflying.set_list_flying_molecules(mylmp=mylmp,region_to_insert=region_to_insert,mymove=mymove,dic_list=self.dic_list,output=mcout)

        #################################################
        # Main Grand canonical MD cycle
        # 1. Equilibrium phase
        # 2. Hybrid NEMD/MC phase for the exchange
        #################################################
        for imc in np.arange(n_MC_s,n_MC_s+n_MC):
            # log file to save all the LAMMPS info
            mylmp.command("log logfile_hneqMDMC.%d"%(imc))
            print ("\n===================\n= EqMD (%d/%d) "%(imc,n_MC),file=mcout)
            mymd.run_MD(mylmp=mylmp,ff=True,md="eqMD",nsteps=n_eqMD,dic_list=self.dic_list,output=mcout,mpiexcept=MPIAbortException)
            if self.dic_list['inputs']['eqMD']['ensemble']=="npt":
                self.dic_list['mc']['volume']=mymisc.get_volume(mylmp=mylmp)
            ###
            # hybrid NEMD/MC for the exchange
            ###
            print ("===================\n= hybrid neqMD/MC (%d/%d)"%(imc,n_MC),file=mcout)
            step_now=mydics.update_current_step(mylmp=mylmp,dic_inputs=self.dic_inputs)
            mymonte.save_old_state(mylmp=mylmp,output=mcout,dic_list=self.dic_list)
            for jmc in np.arange(n_MC_each):
                #mylmp.command("log logfile_hneqMDMC.%d.%d"%(jmc,imc))
                # calculate initial energies in Eq PES (for several sampling)
                mymd.run_MD(mylmp=mylmp,ff=True,md="eqMD",nsteps=0,dic_list=self.dic_list,output=mcout,mpiexcept=MPIAbortException,dump=False,shake=False)
                print ("===================\n= Select a move (%d/%d)"%(jmc,n_MC_each),file=mcout)
                flag_instant=mymonte.prep_MC(mylmp=mylmp,output=mcout,dic_list=self.dic_list,istep=imc)
                # check if an instantaneous MC will be run
                n_neqMD_sel,flag_instant=self.check_if_instant_MC(n_neqMD=n_neqMD,flag_instant=flag_instant,output=mcout)
                print ("===================\n= neqMD to propagate (%d/%d) "%(jmc,n_MC_each),file=mcout)
                mydics.check_current_noneq_dic(output=mcout,dic_inputs=self.dic_inputs)
                flag_neq=mymd.run_MD(mylmp=mylmp, ff=True, md="neqMD", nsteps=n_neqMD_sel, dic_list=self.dic_list, output=mcout, loop=self.dic_list['mc']['loop'], mpiexcept=MPIAbortException, zero_salt_del=zero_salt_del)
                print ("===================\n= Metropolis-type MC (%d/%d) "%(jmc,n_MC_each),file=mcout)
                mymonte.run_MC(mylmp=mylmp,output=mcout,dic_list=self.dic_list,flag_instant=flag_instant,flag_neq=flag_neq)
                mcout.flush()
                # update current step back to before NEMD
                mylmp.command("reset_timestep %d"%step_now)
                #mylmp.command("write_data restart_after_return_%d_of_%d.data pair ij"%(jmc,imc))
                # ungroup
                mylmp.commands_string(mypair.del_groups())
            if imc%save_freq==0:
                mylmp.command("write_data restart_after_MC_%d.data pair ij"%imc)
                # other options; logfile and trajectory
        mylmp.command("write_data restart_end_run.data pair ij")
        print ("===================\n= computation is successfully done", file=mcout) 
        mcout.close()

if __name__=="__main__":
        # my driver for hneMDMC
        mycontrol=LMP_control(rseed=rseed)
        # run hybrid neMD/MC
        mycontrol.run_hneMDMC(
            read_input=True,
            n_preMD=10000,
            n_MC=1,
            n_eqMD=1000,
            n_neqMD=1200, #60000,
            region_to_insert="confined",
            testprint=True,
            filemc="mc_statistics.out")


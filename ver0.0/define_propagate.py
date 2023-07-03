from ctypes import *
import sys,random,math
import numpy as np
import os

#import deepdiff
#import json

#me = 0
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#me = comm.Get_rank()
#nprocs = comm.Get_size()

# import modules
from dictionaries import LMP_dics
mydics=LMP_dics()
from define_lmp_misc import LMP_misc
mymisc=LMP_misc()
from define_pair import LMP_pair
mypair=LMP_pair()
from define_trial_moves import LMP_trial_moves
mymove=LMP_trial_moves()

# Class for hybrid MD/MC
class LMP_propagate():
    def __init__(self):
        self.intro="Define pair potentials"

    def mom_rev_rand(self,dic=None): #,mylmp=None):
        # symmetric two-end momentum reversal
        # need to run random number only once
        momentum=False
        if random.random()<0.5:
            momentum=True
        if not self.dic_mc["gcmc"]:
            momentum=False
        dic['mom_rev']=momentum

    # main driver to run MD (either eq or noneq)
    def run_MD(self,mylmp=None,ff=True,md="eqMD",nsteps=0, dic_list=None,output="test_mc",loop=False,mpiexcept=None,zero_salt_del=False,dump=True,shake=True):
        flag_neq=False
        dic_inputs=dic_list["inputs"]
        self.dic_mc=dic_list["mc"]
        self.dic_peratom=dic_list["peratom"]
        # define a proper input dictionary
        if md not in dic_inputs.keys():
            mymisc.print_out(m="Wrong md keyword", w="ERROR")
        dic=dic_inputs[md]
        # current time step
        if dic['neq']['flag']:
            dic['neq']['current_step_neq'] = mylmp.extract_global("ntimestep",0)
            dic_eq=dic_inputs["eqMD"]
            # momentum reversal random number
            self.mom_rev_rand(dic=dic) #,mylmp=mylmp)
            mylmp.command("run 0 post no # mom rev random number")
            # front-end momentum reversal
            if dic['mom_rev']:
                self.momentum_reversal(mylmp=mylmp, where="front_end")
        # update force field
        if ff:
            mypair.define_force_field(mylmp=mylmp,dic=dic,dic_peratom=self.dic_peratom)
            # implicit wall
            if dic['confinement']['flag']:
                wall=dic['confinement']['wall']
                eps=dic['confinement']['eps']
                sig=dic['confinement']['sig']
                cut=dic['confinement']['cut']
                if wall=="wall/lj1043":
                    spacing=dic['confinement']['spacing']
                else:
                    spacing=1.0 # arbitrary
                mylmp.commands_string(mypair.confinement_wall(on=True,neq=dic['neq']['flag'],wall=wall,eps=eps,sig=sig,cut=cut,spacing=spacing))
        # dynamical rule
        rand=random.randrange(100000000)
        flag,m=self.propagate(mylmp=mylmp,dic=dic,rseed=rand,dump=dump,shake=shake)
        if not flag:
            mymisc.print_out(m=m,w="ERROR")
        ###
        # run MD
        ###
        try:
            mylmp.command("run %d # (neq? %s) MD run"%(nsteps,dic['neq']['flag']))
        except mpiexcept as ae:
            # Single MPI process got killed. This would normally be handled by an MPI abort
            print ("Here is error handling (lmp): this will be re-run with different random number",file=output)
            flag_neq=True
        except Exception as e:
            # All (MPI) processes have reached this error
            print ("Here is error handling: this will be re-run with different random number",file=output)
            flag_neq=True
        # unfix/undump 
        self.unfix(mylmp=mylmp,dic=dic,dump=dump, shake=shake)
        # post processing for non-eq MD - not require fix or dump
        if dic['neq']['flag']:
            # with dictionary of eqMD
            if zero_salt_del:
                mylmp.command("# zero salt case: trial deletion (no coulomb for LJ electrolytes)")
                dic_eq=dic_inputs["eqMD-0"]
            self.post_propagate(mylmp=mylmp,dic=dic,dic_eq=dic_eq,output=output)
        #else:
        #    mylmp.command("group eq_run delete")
        return flag_neq # this is a flag for error

    # unfix MD/undump
    def unfix(self,mylmp=None,dic=None, dump=True, shake=True):
        # undump
        if dic['dump']['flag'] and dump: 
            mylmp.command("undump prod_neq")
        # unshake
        if dic['shake']['flag'] and shake:
            mylmp.commands_string(dic['shake']['unfix'])
        # unfix
        fixid=dic['fixid']
        ensemble=dic['ensemble']
        if ensemble=="nve" or ensemble=="nvt" or ensemble=="npt":
            mylmp.command("unfix    %s" % fixid[0])
        elif ensemble=="langevin": # or ensemble=="rigid/nve":
            mylmp.command("unfix    %s" % fixid[0])
            mylmp.command("unfix    %s" % fixid[1])
        elif ensemble=="rigid/nve":
            mylmp.command("unfix    %s" % fixid[0])
            if len(fixid)>1:
                mylmp.command("unfix    %s" % fixid[1])
        else:
            return False, "Ensemble should be either nve/nvt/npt/langevin"
        ## energy calculation for old config
        #if dic['confinement']['flag']:
        #    mylmp.commands_string(mypair.confinement_wall(on=False))

    # momentum reversal for hybrid MD/MC
    def momentum_reversal(self,mylmp=None,where=None):
        natoms=mylmp.get_natoms()
        v=(3*natoms*c_double)()
        v=mylmp.gather_atoms("v",1,3)
        for i in range(3*natoms):
            v[i] = -v[i]
        mylmp.scatter_atoms("v",1,3,v)    # send it to LAMMPS
        mylmp.command("run 0 post no # momentum reversal applied")            # internal update
        return "Momentum reversed: (%s)"%where

    # dynamical propation in eq or neq
    def propagate(self,mylmp=None,dic=None,rseed=9090,dump=True,shake=True):
        dt=dic['dt']
        fixid=dic['fixid']
        ensemble=dic['ensemble']
        thermo_freq=dic['thermo_freq']
        # set dump
        dd=dic['dump']
        if dd['flag'] and dump:
            mylmp.commands_string("""
dump prod_neq %s custom %d %s %s
dump_modify     prod_neq sort  id append %s
                    """ % (dd['atom'],dd['freq'],dd['file'],dd['style'],dd['append']))
        # Run MD
        if not dic['neq']['flag']:
            # Equilibrium MD
            # shake
            if not shake:
                mylmp.command("# shake is off; should be a dummy run")
            if ensemble!="langevin" and dic['shake']['flag'] and shake:
                mylmp.commands_string(dic['shake']['fix'])
            ensemble=dic['ensemble']
            mylmp.command("timestep %f" % dt)
            mylmp.command("thermo %d" % thermo_freq) # test - will be removed
            if ensemble=="nve":
                mylmp.command("fix  %s   eq_run nve" % fixid[0])
            else:
                temp=dic['temp']
                damp=dic['temp_damp']
                if ensemble=="nvt":
                    mylmp.command("fix  %s   eq_run nvt temp %f %f %f" % (fixid[0], temp, temp, damp))
                elif ensemble=="npt":
                    press=dic['press']
                    damp_p=dic['press_damp']
                    mylmp.command("fix  %s   eq_run npt temp %f %f %f iso %f %f %f" % (fixid[0], temp, temp, damp, press, press, damp_p))
                elif ensemble=="langevin":
                    mylmp.command("fix  %s   eq_run langevin %f %f %f %d" % (fixid[0], temp, temp, damp,rseed))
                    # shake
                    if dic['shake']['flag'] and shake:
                        mylmp.commands_string(dic['shake']['fix'])
                    mylmp.command("fix  %s   eq_run nve" % (fixid[1]))
                else:
                    return False, "Eq ensemble should be either nve/nvt/npt/langevin"
            return True, "Eq MD run ready"
        else:
            # Nonequilibrium MD
            # shake
            if dic['shake']['flag'] and shake:
                mylmp.commands_string(dic['shake']['fix'])
            mylmp.command("timestep %f" % dt)
            # deterministic microcanonical propagation
            if len(fixid)>1:
                mylmp.command("fix  %s  noneq_run_rigid     rigid/nve molecule  # rigid water" % fixid[0])
                mylmp.command("fix  %s  noneq_run           nve                 # nve for non-rigid" % fixid[1])
            else:
                mylmp.command("fix  %s  noneq_run           nve                 # nve for non-rigid" % fixid[0])
            return True, "non-Eq MD run ready"

    # after non-eq MD
    def post_propagate(self,mylmp=None,dic=None,dic_eq=None,output="mc_test",loop=False):
        note="""
        # After propagation
        # insertion: with the flying particle in 3D space (w=0)
        # deletion: with the flying particle at w=w_max, yet it's not removed
        """
        # back-end momentum reversal
        if dic['mom_rev']:    
            self.momentum_reversal(mylmp=mylmp,where="back_end")
        # calculate kinetic energy and update dictionary
        natoms_new=mylmp.extract_global("natoms",0)     # dummy to check
        ke_new=mylmp.extract_compute("thermo_ke",0,0)
        dic['energies']['ke_new']=ke_new
        pe_new=mylmp.extract_compute("thermo_pe",0,0)   # dummy to check
        print ("new state after neqMD (natom=%d,pe_new=%f,ke_new=%f)"%(natoms_new,pe_new,ke_new),file=output)
        ###
        # return to equilibrium PES with non-interacting flying molecules
        if dic_eq['confinement']['flag']:
            wall=dic_eq['confinement']['wall']
            eps=dic_eq['confinement']['eps']
            sig=dic_eq['confinement']['sig']
            cut=dic_eq['confinement']['cut']
            if wall=="wall/lj1043":
                spacing=dic_eq['confinement']['spacing']
            else:
                spacing=1.0 # arbitrary
            # eq PES - confinement
            if not self.dic_mc['select_insert']: # non-interacting flying molecule
                group="noneq_run"
            else:
                group="all"
            mylmp.commands_string(mypair.confinement_wall(on=True,group=group,neq=dic_eq['neq']['flag'],wall=wall,eps=eps,sig=sig,cut=cut,spacing=spacing))
        #    mylmp.command("run 0 post no # for energy after confinement")
        # return to equilibrium PES
        mypair.define_force_field(mylmp=mylmp,dic=dic_eq,dic_peratom=self.dic_peratom)
        if not self.dic_mc['select_insert']:
            # non-interacting flying molecules
            ljcut=dic_eq['cutoff']['lj']    # equilibrium cutoff
            m=mymove.delete_molecule(mylmp=mylmp, where="trial_delete_after_neMD", ljcut=ljcut, ff=dic_eq['ff'])
            print (m,file=output)
        else:
            mylmp.command("run 0 post no # update forcefield; return to Eq")
        ###
        # calculate potential energy and update dictionary
        pe_new = mylmp.extract_compute("thermo_pe",0,0) 
        # compute intramolecular energies
        bond_new,angle_new=mymisc.compute_water_intramolecular(mylmp=mylmp)
        dic['energies']['bond_w_new']=bond_new
        dic['energies']['angle_w_new']=angle_new
        pe_new-=bond_new+angle_new
        dic['energies']['pe_new']=pe_new                        # after substracting artifacial intramolecular potential
        ke_new = mylmp.extract_compute("thermo_ke",0,0) # dummy to check
        if dic_eq['confinement']['flag']:
            mylmp.commands_string(mypair.confinement_wall(on=False))
        print ("new state after back to eqPES (natom=%d,pe_new=%f,ke_new=%f)"%(natoms_new,pe_new,ke_new),file=output)
        print ("energy chage in water intermolecular energies: bond=(%f -> %f) / angle=(%f -> %f)"%(dic['energies']['bond_w_old'],bond_new,dic['energies']['angle_w_old'],angle_new),file=output)
        mylmp.command("# done post_propagate after NEMD")

    def instant_delete(self,mylmp=None,dic=None):
        self.delete_molecule(mylmp=mylmp,where="instant_delete")
        # calculate kinetic energy
        ke_new = mylmp.extract_compute("thermo_ke",0,0)
        dic['energies']['ke_new']=ke_new

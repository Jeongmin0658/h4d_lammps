from ctypes import *
import sys,random,math
import numpy as np
import os

# Class for hybrid MD/MC
class LMP_misc():
    def __init__(self):
        self.intro="Define LAMMPS misc"
        self.mass_oxygen=15.999400
        self.mass_hydrogen=1.007940
        self.mass_diff=(self.mass_oxygen-self.mass_hydrogen)
        self.mass_uniform=(self.mass_oxygen+2*self.mass_hydrogen)/3.

        # units
        self.kJ_to_kcal=1.0/4.184
        self.kT_cal=0.593       # thermal energy in kcal/mol
        self.Avogadro=6.022

    def print_out(self,mylmp=None,m="Default message",w="hMDMC"):
        mylmp.command("print '%s ::: %s'"%(w,m))
        if w=="ERROR":
            mylmp.command("print '%s ::: Terminated unsuccefully'" % w)
            exit()

    def calc_thermo(self,mylmp=None):
        mylmp.command("compute thermo_ke all ke")
        #mylmp.command("compute thermo_pe all pe")

    # compute intramolecular energies
    def compute_water_intramolecular(self,mylmp=None):
        #mylmp.command("run 0 post no # to compute water intermolecular energies")
        bond_old=mylmp.get_thermo("ebond")
        angle_old=mylmp.get_thermo("eangle")
        return bond_old,angle_old

#    def define_groups(self,mylmp=None,group=None):
#        # define groups of each species
#        #mylmp.command("group oxygens    type 1")
#        #mylmp.command("group hydrogens  type 2")
#        #mylmp.command("group waters     type 1 2")
#        #mylmp.command("group cations    type 3")
#        #mylmp.command("group anions     type 4")
#        #mylmp.command("group flying_waters      type 8 9")
#        #mylmp.command("group flying_salts       type 10 11")
#        mylmp.command("group all_waters type 1 2 8 9")
#        mylmp.commands_string(group)

    ### begin - per-atom properties at current state
    # allocate per-atom properties: position, velocity, force, images
    def allocate_xvf(self,mylmp=None):
        natoms=mylmp.get_natoms()
        x=(3*natoms*c_double)()
        v=(3*natoms*c_double)()
        f=(3*natoms*c_double)()
        q=(natoms*c_double)()
        atype=(natoms*c_int)()
        images=(3*natoms*c_int)()
        return x,v,f,q,atype,images

    def get_volume(self,mylmp=None):
        boxlo,boxhi,xy,yz,xz,periodicity,box_change=mylmp.extract_box()
        self.volume=1.
        for j in range(3):
            self.volume*=boxhi[j]-boxlo[j]
        return self.volume

    def gather_atoms_xvf(self,mylmp=None):
        x,v,q,atype,f,images=self.allocate_xvf(mylmp=mylmp)
        x = mylmp.gather_atoms("x",1,3) # wrapped
        v = mylmp.gather_atoms("v",1,3)
        f = mylmp.gather_atoms("f",1,3)
        q = mylmp.gather_atoms("q",1,1)
        #mylmp.command("# after gather x,v,f,q")
        atype = mylmp.gather_atoms("type",0,1)
        images = mylmp.gather_atoms("image",0,3)
        return x,v,f,q,atype,images
    ### end - per-atom properties at current state
#
#    def init_peratom(self,mylmp=None,dic_peratom=None):
#        x,v,f,image=self.allocate_xvf(mylmp=mylmp)
#        for item in ["new","old","current"]:
#            dic_peratom[item]['x']=x
#            dic_peratom[item]['v']=v
#            dic_peratom[item]['f']=f
#            dic_peratom[item]['image']=image
#
#    def update_peratom(self,mylmp=None,dic_peratom=None,status="old"):
#        # update dictionary
#        x_old,v_old,f_old,image_old=self.gather_atoms_xvf(mylmp=mylmp)
#        dic_peratom[status]['x']=x_old
#        dic_peratom[status]['v']=v_old
#        dic_peratom[status]['f']=f_old
#        dic_peratom[status]['image']=image_old

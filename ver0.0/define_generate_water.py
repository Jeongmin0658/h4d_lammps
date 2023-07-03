from ctypes import *
import sys,random,math
import numpy as np
import os

# import modules
from define_flying_water_init import LMP_flying_init
myflyinit=LMP_flying_init()

# Class for hybrid MD/MC
class LMP_generate_water():
    def __init__(self):
        self.intro="Define water generation"

        # SPC/E water
        # angle in radian
        self.a_oh=0.5*(180.-109.47)/180.*math.pi
        self.d_oh=1.  # angstrom
        # two hydrogens
        ha=math.cos(self.a_oh)*self.d_oh
        hc=math.sin(self.a_oh)*self.d_oh
        # residual vector
        x_rel=[0.,0.,0.]
        x_rel.extend([ha,0.,hc])
        x_rel.extend([-ha,0.,hc])
        self.x_rel=x_rel

    # generate a single water molecule
    def generate_water(self,mylmp=None,dic_inputs=None,mass=[16,1,1],boxlo=[-1.,-1.,-1.],boxhi=[1.,1.,1.],test_structure=False):
        ###
        # configuration space
        ###
        # oxygen - center of mass
        pos_o=[]
        for j in np.arange(3):
            pos_o.append(random.uniform(boxlo[j],boxhi[j]))
        # rotate the vector
        x_new=myflyinit.rot_water(self.x_rel)    # new orientation
        mylmp.command("# x_rel_new: %f %f %f %f %f %f %f %f %f"%tuple(x_new))
        # check water structure - checked already
        if test_structure:
            l=[]
            for xi in [x_new[3:6],x_new[6:9]]:
                n=np.dot(xi,xi)
                n=math.sqrt(n)
                l.append(n)
                mylmp.command("# x_rel_new_norm: %f"%n)
            hh=np.subtract(x_new[3:6],x_new[6:9])
            l.append(math.sqrt(np.dot(hh,hh)))
            a=np.arccos((l[0]*l[0]+l[1]*l[1]-l[2]*l[2])/(2.*l[0]*l[1]))*180./math.pi
            mylmp.command("# x_rel_new_angle: %f"%a)
        # space fame
        oxy=pos_o[:]
        oxy.extend(pos_o)
        oxy.extend(pos_o)
        x_new=np.add(x_new,oxy)            # back to space frame
        mylmp.command("# x_new_before_pbc: %f %f %f %f %f %f %f %f %f"%tuple(x_new))
        # pbc check for hydrogen
        box=[boxhi[j]-boxlo[j] for j in np.arange(3)]
        mylmp.command("# box size: %f %f %f"%tuple(box))
        # image of hydrogens
        imageh=[[round(x/b) for x,b in zip(x_new[3:6],box)],[round(x/b) for x,b in zip(x_new[6:9],box)]]
        # wrapped xyz of hydrogens
        pos_h1=[x-round(x/b)*b for x,b in zip(x_new[3:6],box)]
        pos_h2=[x-round(x/b)*b for x,b in zip(x_new[6:9],box)]
        ###
        # momentum space
        ###
        # thermal energy
        trans=dic_inputs['thermal']/10000. # g/mol*(A/s)^2
        rot=dic_inputs['thermal']/10000. # g/mol*(A/s)^2
        # get initialized velocity
        v_new=myflyinit.random_atom_vel_rigid(mass=mass,x=x_new,trans=trans,rot=rot)
        mylmp.command("# v_new: %f %f %f %f %f %f %f %f %f"%tuple(v_new))
        ## update velocity array
        ## need to be corrected
        #natoms=mylmp.get_natoms()
        #v_f=(3*natoms*c_double)()
        #v_f=mylmp.gather_atoms("v",1,3)
        #for j in np.arange(9):
        #    v_f[-9+j]=v_new[j]
        #mylmp.scatter_atoms("v",1,3,v_f)            # send it to LAMMPS
        #mylmp.command("run 0 post no # new flying water velocity updated")
        return pos_o,[pos_h1,pos_h2],imageh,v_new


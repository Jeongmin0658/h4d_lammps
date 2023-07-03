from ctypes import *
import sys,random,math
import numpy as np
import os
import timeit
#from scipy import integrate

## Dictionary as input files
#from dictionaries import LMP_dics
#mydics=LMP_dics()

##rank = 0
#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()

# Class for hybrid MD/MC
class LMP_define_select_move():
    def __init__(self):
        self.intro="Define a set of water to move while a switching"

    def func_dist_calc(self, d=0., cut=7., w=0.1):
        # d: distance from a flying ion
        return 0.5*(1.-np.tanh((d-cut)/w))

    def select_criteria(self, rp=0., rm=0.):
        # rp: distance from flying cation
        # rm: distance from flying anion
        return 1.-(1.-rp)*(1.-rm)

    def find_free_noneq(self, mylmp=None, x=None, box=None, ci=None, dic_mc=None, output="mc_out", cut=7, width=0.1):
        # flying ions
        if dic_mc["select_salt"]:
            cation=dic_mc['selected']['salt']['cation']-1
            anion=dic_mc['selected']['salt']['anion']-1
            xc=x[3*cation:3*cation+3]           # cation position
            xa=x[3*anion:3*anion+3]             # anion position
            # getting particle indices to relax
            sel_list=[]
            for c in ci:
                xq=x[3*c:3*c+3] # position vector
                sel=self.get_factor_from_dist(xq,xc,xa,box,cut=cut, width=width)
                sel_list.append(sel)
            return np.array(sel_list),ci # list of particles to relax and its factor
        elif dic_mc["select_water"]:
            water=dic_mc['selected']['water']-1
            xw=x[3*water:3*water+3]           # cation position
            mylmp.command("# init selection center (water): %s"%str(xw))
            # getting particle indices to relax
            sel_list=[]
            sel_w=[]
            for c in dic_mc['prop_to_save']['water_list']:
                if c-1!=water:
                    xq=x[3*c-3:3*c] # position vector
                    sel=self.get_factor_from_dist_single(xq,xw,box,cut=cut, width=width)
                    sel_list.append(sel)
                    sel_w.append(c-1)
            return np.array(sel_list),np.array(sel_w) # list of particles to relax and its factor
        else:
            exit("LJ is not ready for find_free_noneq")

    def select_process(self, ci, sel_list):
        # draw a random number and select
        free_move=[]
        for c,sel in zip(ci,sel_list):
            rd=random.random()
            if sel>=rd: # select free moving particles
                free_move.append(c+1)
        return np.array(free_move)

    # for single water or LJ
    def get_factor_from_dist_single(self,xq,xc,box,cut=7, width=0.1):
        xqc=np.subtract(xq,xc)
        xqc=np.subtract(xqc, box*np.rint(np.divide(xqc,box)))
        lc=np.sqrt(xqc.dot(xqc))
        sel=self.func_dist_calc(d=lc, cut=cut, w=width) # tanh
        return sel

    # for ion pair
    def get_factor_from_dist(self,xq,xc,xa,box,cut=7, width=0.1):
        #print ("array vector size:",xq,xc,xa,"rank=",rank)
        xqc=np.subtract(xq,xc)
        xqa=np.subtract(xq,xa)
        # pbc
        xqc=np.subtract(xqc, box*np.rint(np.divide(xqc,box)))
        xqa=np.subtract(xqa, box*np.rint(np.divide(xqa,box)))
        # get distances from either cation or anion
        lc=np.sqrt(xqc.dot(xqc))
        la=np.sqrt(xqa.dot(xqa))
        # activation function; hyperbolic tangent
        sc=self.func_dist_calc(d=lc, cut=cut, w=width)
        sa=self.func_dist_calc(d=la, cut=cut, w=width)
        # final check: highy likely overlap 
        sel=self.select_criteria(rp=sc, rm=sa)
        return sel

    # factor calculation for acceptance probability
    def get_factor_selection(self, mylmp=None, dic_mc=None, output="mc_out", cut=7, width=0.1):
        natoms=mylmp.get_natoms()
        x=(3*natoms*c_double)()
        x=mylmp.gather_atoms("x",1,3) # wrapped position
        boxl=mylmp.extract_box()
        box=np.subtract(boxl[1],boxl[0])    # box size
        if dic_mc['select_salt']:
            # flying ions
            cation=dic_mc['selected']['salt']['cation']-1
            anion=dic_mc['selected']['salt']['anion']-1
            xc=x[3*cation:3*cation+3]           # cation position
            xa=x[3*anion:3*anion+3]             # anion position
            # initialize factor
            factor=1.
            mylmp.command("run 0 pre no post no # getting a factor relevant for selection (after init; factor=%f)"%factor)
            for c in self.free_list:
                xq=x[3*c-3:3*c] # c = c-1; position vector
                sel=self.get_factor_from_dist(xq,xc,xa,box, cut=cut, width=width)
                factor*=sel
        elif dic_mc['select_water']:
            # flying water
            water=dic_mc['selected']['water']-1
            xw=x[3*water:3*water+3]           # cation position
            # initialize factor
            factor=1.
            mylmp.command("run 0 pre no post no # getting a factor relevant for selection (after init; factor=%f; center=%s)"%(factor,str(xw)))
            #print ("# again free_list selection (water, %f/%f):"%(cut,width),len(self.free_list),self.free_list,file=output)
            for c in self.free_list:
                xq=x[3*c-3:3*c] # c = c-1; position vector
                sel=self.get_factor_from_dist_single(xq,xw,box, cut=cut, width=width)
                factor*=sel
        else:
            mylmp.command("# LJ is picked for selection:::under construction")
            print ("Not ready now:::get_factor_selection")
            exit()
        mylmp.command("run 0 pre no post no # getting a factor relevant for selection (factor=%f;%d)"%(factor,len(self.free_list)))
        return factor

    def select_free_to_move(self, mylmp=None, dic_mc=None, output="mc_out", cut=7, width=0.1):
        # now get the list of free atoms
        natoms=mylmp.get_natoms()
        x=(3*natoms*c_double)()
        x=mylmp.gather_atoms("x",1,3) # wrapped position
        boxl=mylmp.extract_box()
        box=np.subtract(boxl[1],boxl[0])    # box size
        # list of all the atoms
        cip = range(natoms)
        # get a list of prob. to select
        free_factor,cip=self.find_free_noneq(mylmp=mylmp, x=x, box=box, ci=cip, dic_mc=dic_mc, output="mc_out", cut=cut, width=width)
        # selection process
        free=self.select_process(cip,free_factor)
        self.free_list=free   # for factor calculation later both before and after NEMD
        #print ("List of particle to relax (tot=%d)"%natoms,len(self.free_list),self.free_list, file=output)
        #print ("List of prob to select",len(free_factor),natoms,free_factor, file=output)
        mylmp.command("run 0 pre yes post no # after selecting free particles")
        myfree="group noneq_free id"
        for i in free:
            myfree+=" %d"%(i)
        mylmp.command("run 0 pre no post no # after selecting free particles")
        mylmp.command(myfree)
        mylmp.command("group noneq_free include molecule")
        mylmp.command("group noneq_free_water intersect   noneq_free water")
        mylmp.command("group noneq_free_other  subtract   noneq_free water")
        mylmp.command("run 0 pre yes post no # after selecting free particles")

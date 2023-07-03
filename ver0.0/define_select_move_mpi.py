from ctypes import *
import sys,random,math
import numpy as np
import os
import timeit
#from scipy import integrate

## Dictionary as input files
#from dictionaries import LMP_dics
#mydics=LMP_dics()

#rank = 0
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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
            cation=dic_mc['selected']['salt']['cation']-1
            anion=dic_mc['selected']['salt']['anion']-1
            xc=x[3*cation:3*cation+3]           # cation position
            xa=x[3*anion:3*anion+3]             # anion position
            # getting particle indices to relax
            free_move=[]
            free_factor=[]
            for c in ci:
                xq=x[3*c:3*c+3]
                #print ("array vector size:",xq,xc,xa,'atom index',c,"rank=",rank)
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
                # draw random number and check
                rd=random.random()
                if sel>=rd: # select free moving particles
                    #print (c+1,sel,rd,s,'check',file=output)
                    free_move.append(c+1)
                    free_factor.append(sel)
            #print ("# Free_move array in proc=%d:"%rank, len(free_move), free_move)
            return np.array(free_move),np.array(free_factor)

    def get_factor_selection(self, mylmp=None): #, dic_mc=None): #, output="mc_out", cut=7, width=0.1):
        if rank == 0:
            #print ("# Total free atoms (get_factor_selection)",len(self.free_list),self.free_list)
            factor=1.
            for s in range(size):
                for f in self.free_list[s]:
                    factor*=f
        else:
            factor=None
        factor=comm.bcast(factor, root=0)
        mylmp.command("run 0 post no # getting a factor relevant for selection (factor=%f)"%factor)
        return factor

    def select_free_to_move(self, mylmp=None, dic_mc=None, output="mc_out", cut=7, width=0.1):
        # now get the variable "selca" and get the list of free atoms
        natoms=mylmp.get_natoms()
        x=(3*natoms*c_double)()
        x=mylmp.gather_atoms("x",1,3) # wrapped position
        boxl=mylmp.extract_box()
        box=np.subtract(boxl[1],boxl[0])    # box size
        # master process; 
        # ref: https://www.kth.se/blogs/pdc/2019/08/parallel-programming-in-python-mpi4py-part-1/
        # distribute x into each process
        if rank == 0:
            # determine the size of each sub-task
            nprocs=size
            ave, res = divmod(natoms, nprocs)
            counts = [ave + 1 if p < res else ave for p in range(nprocs)]
            t = sum(counts)
            print (ave,res,counts,'divmod and total=natoms',t,natoms)
        
            # determine the starting and ending indices of each sub-task
            starts = [sum(counts[:p]) for p in range(nprocs)]
            ends = [sum(counts[:p+1]) for p in range(nprocs)]
            print ("starts  :",starts)
            print ("ends    :",ends)
        
            # converts data into a list of arrays 
            #xip = [xi[starts[p]:ends[p]] for p in range(nprocs)]
            cip = [range(starts[p],ends[p]) for p in range(nprocs)]
            #cxip = [[x,c] for x,c in zip(xip,cip)]
        else:
            #cxip = None
            cip = None
        cip = comm.scatter(cip, root=0)
        print('#Process {} has data:'.format(rank), len(cip), "first item", cip[0], file=output)
        free,free_factor=self.find_free_noneq(mylmp=mylmp, x=x, box=box, ci=cip, dic_mc=dic_mc, output="mc_out", cut=7, width=0.1)
        # broadcast
        free=comm.gather(free, root=0)
        free_factor=comm.gather(free_factor, root=0)
        self.free_list=free_factor   # for factor calculation later both before and after NEMD
        if rank==0:
            #print ("# Total free atoms",len(free),free)
            myfree="group noneq_free id"
            for s in range(size):
                for i in free[s]:
                    myfree+=" %d"%(i)
        else:
            myfree=None
        myfree=comm.bcast(myfree, root=0)
        mylmp.command(myfree)
        mylmp.command("group noneq_free include molecule")
        mylmp.command("group noneq_free_water intersect   noneq_free water")
        mylmp.command("group noneq_free_other  subtract   noneq_free water")
        mylmp.command("run 0 post no # after selecting free particles")
        #exit()
        #mylmp.command("quit")
        #end = timeit.timeit()

### Not use below
        # flying ions
#        cation=dic_mc['selected']['salt']['cation']
#        anion=dic_mc['selected']['salt']['anion']
#        seed=random.randrange(1,1000000) # random seed
#        mylmp.commands_string("""
####
## selection process for free moving particles
####
## reference point: 3D positions of flying cation/anion
## hyperbolic tangent with cut={cut} and width={width}
#### from cation
#variable cx     atom x-x[{cation}]
#variable cy     atom y-y[{cation}]
#variable cz     atom z-z[{cation}]
## pbc
#variable cxp    atom v_cx-round(v_cx/lx)*lx
#variable cyp    atom v_cy-round(v_cy/ly)*ly
#variable czp    atom v_cz-round(v_cz/lz)*lz
## dist from cation
#variable distc  atom sqrt(v_cxp*v_cxp+v_cyp*v_cyp+v_czp*v_czp)
#variable distce atom exp(2.0*(v_distc-{cut})/{width})   # arg=(r-R)/w
#variable selc   atom 0.5*(1.0-(v_distce-1.0)/(v_distce+1.0)) # hyperbolic
#### from anion
#variable ax     atom x-x[{anion}]
#variable ay     atom y-y[{anion}]
#variable az     atom z-z[{anion}]
## pbc
#variable axp    atom v_ax-round(v_ax/lx)*lx
#variable ayp    atom v_ay-round(v_ay/ly)*ly
#variable azp    atom v_az-round(v_az/lz)*lz
## dist from anion
#variable dista  atom sqrt(v_axp*v_axp+v_ayp*v_ayp+v_azp*v_azp)
#variable distae atom exp(2.0*(v_dista-{cut})/{width})   # arg=(r-R)/w
#variable sela   atom 0.5*(1.0-(v_distae-1)/(v_distae+1.0)) # hyperbolic
#### selection process from two selection functions
## variable selca should be either 0.0 (false) or 1.0 (true)
#variable selca  atom 1.0-(1.0-v_selc)*(1.0-v_sela) #>= random(0,1,{seed})
#run 0 post no # getting selection variable
#        """.format(cation=cation, anion=anion, seed=seed, cut=cut, width=width)
#        )

        #mylmp.command("run 0 post no # just before selection")
        #if rank == 0:
        #    myfree="group noneq_free id"
        #    count,mv=0,0
        #    for i,sel in enumerate(selca):
        #        rd=random.random()
        #        if sel>=rd: # select free moving particles
        #            #print (i+1,sel,rd,s,'check',file=output)
        #            myfree+=" %d"%(i+1)
        #            mv+=1
        #        count+=1
        #    myfree+=" # (%d) total number="%mv+str(count)+"/"+str(natoms)
        #    print (myfree,'inside rank 0',file=output)
        #else:
        #    myfree=None
        ###

    def remove_vars(self, mylmp=None):
        mylmp.commands_string("""
### remove variables used in selection processes
variable cx     delete
variable cy     delete
variable cz     delete
variable cxp    delete
variable cyp    delete
variable czp    delete
variable distc  delete
variable distce delete
variable selc   delete
variable ax     delete
variable ay     delete
variable az     delete
variable axp    delete
variable ayp    delete
variable azp    delete
variable dista  delete
variable distae delete
variable sela   delete
variable selca  delete
###
        """)

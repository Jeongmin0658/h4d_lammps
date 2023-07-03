from ctypes import *
import sys,random,math
import numpy as np
import os
#from scipy import integrate

## Dictionary as input files
#from dictionaries import LMP_dics
#mydics=LMP_dics()

# Class for hybrid MD/MC
class LMP_define_bias_salt():
    def __init__(self):
        self.intro="Define bias function for ion-pair exchange"
        # Salt pair
        # bias function: either bimodal or gaussian

    ### bias function without excluded volume
    def prob_bias_select_instant(self,box,alpha_bias,rcut,x_center=1.,bias_f="bimodal"):
        # need to be fix; 17/May/2022
        if bias_f=="gaussian":
            prob=math.erf(math.sqrt(alpha_bias)*rcut)/math.erf(math.sqrt(alpha_bias)*box*0.5)
        elif bias_f=="bimodal":
            norm=math.erf(math.sqrt(alpha_bias)*(x_center+box*0.5))-math.erf(math.sqrt(alpha_bias)*(x_center-box*0.5))
            norm_cut=math.erf(math.sqrt(alpha_bias)*(x_center+rcut))-math.erf(math.sqrt(alpha_bias)*(x_center-rcut))
            prob=norm_cut/norm
        elif bias_f=="uniform":
            norm=box
            norm_cut=rcut*2.
            prob=norm_cut/norm
        else:
            print ("ERROR:::Not proper bias function is selected")
            exit()
        return np.power(prob,3.)  # xyz

    def func_dist_bias(self,x=0.,box=1.,bias_f="bimodal",alpha_bias=1.,x_center=1.):
        # target dist
        if bias_f=="gaussian":
            norm=math.sqrt(alpha_bias/math.pi)/math.erf(math.sqrt(alpha_bias)*box*0.5)
            prob=norm*math.exp(-alpha_bias*x*x)
        elif bias_f=="bimodal":
            norm=math.sqrt(alpha_bias/math.pi)/(math.erf(math.sqrt(alpha_bias)*(x_center+box*0.5))-math.erf(math.sqrt(alpha_bias)*(x_center-box*0.5)))
            prob=norm*(math.exp(-alpha_bias*np.power(x-x_center,2.))+math.exp(-alpha_bias*np.power(x+x_center,2.)))
        elif bias_f=="uniform":
            norm=1./box # test
            prob=norm
        else:
            print ("ERROR:::Not proper bias function is selected")
            exit()
        #ratio=1. # test
        return prob,norm

    def func_dist_uniform(self,x=0.,box=1.,bias_f="bimodal",alpha_bias=1.,x_center=1.):
        # target dist
        prob,const=self.func_dist_bias(x=x,box=box,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
        # the ratio
        if bias_f=="bimodal":  # max probability; it is 1 for Gaussian
            const*=(1.+math.exp(-4.*alpha_bias*x_center))
        ratio=prob/const # the ratio using random number
        #ratio=1. # test
        if bias_f=="uniform" and ratio != 1.:
            exit("ERROR::: uniform bias function should be ratio=1")
        rd=-random.random()+1.      # random number
        if rd <= ratio:
            return True,prob
        return False,prob
    ###
        
    def bias_after_NEMD(self,mylmp=None,x=None,forward=True,dic_mc=None):
        ## flying ions
        cation=dic_mc['selected']['salt']['cation']
        anion=dic_mc['selected']['salt']['anion']
        anions_list=dic_mc['prop_to_save']['anions_list']
        # bias functions
        bias_f=dic_mc['bias']['bias_f']
        alpha_bias=dic_mc['bias']['alpha']
        x_center=dic_mc['bias']['center']
        mylmp.command("# selected salt %d/%d"%(cation,anion))
        mylmp.command("#The center of bimodal dist. on a positive side=%f"%x_center)
        # box for pbc
        boxlo,boxhi,xy,yz,xz,periodicity,box_change=mylmp.extract_box()
        box=[boxhi[j]-boxlo[j] for j in np.arange(3)]
        if forward: # after insertion
            atom_c=3*(cation-1)
            pos_c=[x[atom_c],x[atom_c+1],x[atom_c+2]]
            # bias distribution of interatomic distance - 
            pw=0.
            for ndelete_a in anions_list:
                atom_a=3*(ndelete_a-1)
                tri_a=[x[atom_a],x[atom_a+1],x[atom_a+2]]
                p=self.bias_salt_del_choice(mylmp=mylmp,salt=[pos_c,tri_a],bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                pw+=p
            # chosen flying anion
            atom_a=3*(anion-1)
            pos_a=[x[atom_a],x[atom_a+1],x[atom_a+2]]
            # wrapped interionic vector
            d=np.subtract(pos_a,pos_c)
            dvec=[x-round(x/b)*b for x,b in zip(d,box)]
            # biased insertion of anion
            prob=1.
            for d,b in zip(dvec,box):
                # bias function
                flag,p=self.func_dist_uniform(x=d,box=b,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                #p=1. # test
                prob*=p
            pw+=prob    # this is not in the array
            # update the factor 
            dic_mc["bias"]["factor_del"]=pw/prob # N+1 if equally probable
        else:   # after deletion
            # calculate prob of accepting this interionic distance to insert
            atom_c=3*(cation-1)
            pos_c=[x[atom_c],x[atom_c+1],x[atom_c+2]]
            atom_a=3*(anion-1)
            pos_a=[x[atom_a],x[atom_a+1],x[atom_a+2]]
            # wrapped interionic vector
            d=np.subtract(pos_a,pos_c)
            dvec=[x-round(x/b)*b for x,b in zip(d,box)]
            # bias function
            prob=1.
            for d,b in zip(dvec,box):
                p,n=self.func_dist_bias(x=d,box=b,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                prob*=p
            # update the factor
            dic_mc['bias']["factor_ins"]=1./prob  # Volume if random insertion

    def bias_salt_del_choice(self,mylmp=None,salt=[[0.,0.,0.],[0.,0.,0.]],bias_f="bimodal",alpha_bias=1.0,x_center=1.0):
        self.alpha_bias=alpha_bias
        # wrapped interionic vector
        boxlo,boxhi,xy,yz,xz,periodicity,box_change=mylmp.extract_box()
        box=[boxhi[j]-boxlo[j] for j in np.arange(3)]
        # calculate prob of accepting this interionic distance of deleting pair
        d=np.subtract(salt[1],salt[0])
        # wrapped interionic vector
        d=[x-round(x/b)*b for x,b in zip(d,box)]
        # bias function
        prob=1.
        for x,b in zip(d,box):
            p,n=self.func_dist_bias(x=x,box=b,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
            prob*=p
        #prob=1. # test
        return prob # probability of accepting this interionic distance

    # INSERTION
    def generate_salt(self,mylmp=None,boxlo=[-1.,-1.,-1.],boxhi=[1.,1.,1.],bias_f="bimodal",alpha_bias=0.01,x_center=1.,rcut=1.):
        pos_c,pos_a,prob,image_a,f_instant=self.bias_salt_choice(mylmp=mylmp,salt=[[],[]],select=False,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center,rcut=rcut)
        return pos_c,pos_a,prob,image_a,f_instant

    def bias_salt_choice(self,mylmp=None,salt=[[0.,0.,0.],[0.,0.,0.]],select=False,bias_f="bimodal",alpha_bias=1.,x_center=1.,rcut=1.):
        # wrapped interionic vector
        boxlo,boxhi,xy,yz,xz,periodicity,box_change=mylmp.extract_box()
        box=[boxhi[j]-boxlo[j] for j in np.arange(3)]
        if select:
        # calculate prob of accepting this interionic distance of deleting pair
            d=np.subtract(salt[1],salt[0])
            # wrapped interionic vector
            d=[x-round(x/b)*b for x,b in zip(d,box)]
            # bias function
            prob=1.
            for x,b in zip(d,box):
                p,n=self.func_dist_bias(x=x,box=b,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                prob*=p
            rd=-random.random()+1.
            #prob=1. # test
            if rd<=prob:
                return True,prob
            return False,prob
        else:
        # generate a pair of cation/anion
            # random insertion of cation
            pos_c=[]
            for j in np.arange(3):
                pos_c.append(random.uniform(boxlo[j],boxhi[j]))
            # biased insertion of anion
            pos_a=[]
            prob=1.
            f_instant=True
            while len(pos_a)<3:
                # interionic vec
                x=random.uniform(boxlo[len(pos_a)],boxhi[len(pos_a)])
                # bias function
                flag,p=self.func_dist_uniform(x=x,box=box[len(pos_a)],bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                if flag:
                    tri=pos_c[len(pos_a)]+x # cation position + interionic vec
                    pos_a.append(tri)
                    prob*=p                 # probability of biased insertion
                    # check if this is in the excluded volume of a cubic box
                    if abs(x)<rcut:
                        if f_instant:
                            f_instant=True
                    else:
                        f_instant=False
            # pbc applied for anion
            image_a=[round(x/b) for x,b in zip(pos_a,box)]
            # wrapped xyz of anion 
            pos_a=[x-round(x/b)*b for x,b in zip(pos_a,box)]
        return pos_c,pos_a,prob,image_a,f_instant

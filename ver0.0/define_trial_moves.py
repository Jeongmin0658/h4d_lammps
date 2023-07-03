from ctypes import *
import sys,random,math
import numpy as np
import os

# Dictionary as input files
from dictionaries import LMP_dics
mydics=LMP_dics()

# import modules
from define_pair import LMP_pair
mypair=LMP_pair()
from define_flying_water_init import LMP_flying_init
myflyinit=LMP_flying_init()

# Need to import bias_salt/water
from define_generate_water import LMP_generate_water
mywater=LMP_generate_water()
from define_bias_salt import LMP_define_bias_salt
mysalt=LMP_define_bias_salt()

# Class for hybrid MD/MC
class LMP_trial_moves():
    def __init__(self):
        self.intro="Define trial moves"

    # rejected move: return to old configuration
    def return_position(self,mylmp=None,dic_peratom=None):
        mylmp.scatter_atoms("q",1,1,dic_peratom['old']['q'])                # send it to LAMMPS
        mylmp.scatter_atoms("x",1,3,dic_peratom['old']['x'])                # send it to LAMMPS
        mylmp.scatter_atoms("v",1,3,dic_peratom['old']['v'])                # send it to LAMMPS
        mylmp.scatter_atoms("f",1,3,dic_peratom['old']['f'])                # send it to LAMMPS
        mylmp.scatter_atoms("image",0,3,dic_peratom['old']['image'])        # send it to LAMMPS
        mylmp.command("run 0 post no # return_position back to old")        # internal update
        return "Return to the old state"

    # accept trial insertion
    def accept_insert(self,mylmp=None,dic_selected=None, where="accept_insert"):
        # types: from noneq to eq (to make flying molecule indistinguishable)
        for i,j in zip(dic_selected['type_neq_l'],dic_selected['type_eq_l']):
            mylmp.command("set type %d type %d" % (i,j))
        mylmp.command("run 0 post no # lj type change (%s)"%where)            # internal update
        return "Accept inserting flying particles (now indistinguishable): (%s)"%(where)

    def type_change_molecule(self,mylmp=None,dic_mc=None,dic_selected=None,where="delete_trial"):
        for i,j in zip(dic_mc['selected']['list'],dic_selected['full_list']['type_neq_l']):
            mylmp.command("set atom %d type %d" % (i,j))
        #mylmp.command("run 0 post no ")            # internal update
        ## charge - no need
        #for i,j in zip(dic_selected['type_neq_l'],dic_selected['charge']):
        #    mylmp.command("set type %d charge %f" % (i,j))
        return"Type change flying particles: (%s)"%(where)

    def charge_return(self,mylmp=None,dic_mc=None,dic_selected=None,where="reject_delete"):
        for i,j in zip(dic_mc['selected']['list'],dic_selected['full_list']['charge']):
            mylmp.command("set atom %d charge %f" % (i,j))
        mylmp.command("run 0 post no # charge return (%s)"%where)            # internal update
        return"Charge return flying particles: (%s)"%(where)

    # velocity generation for flying ions
    def init_vel_single(self,mylmp=None,dic=None,where=None,temp=300.):
        rseed1=random.randint(1,10000000)
        mylmp.command("velocity widom_test create %f %d dist gaussian mom no rot no units box" % (temp,rseed1))    # send it to LAMMPS
        mylmp.command("run 0 post no # temp=%f velocity generated"%temp)

    def init_vel_single_lj(self,mylmp=None,dic=None,where=None,temp=300.,atom_index=0, mass=1., unit_lj=True):
        v = mylmp.gather_atoms("v",1,3)
        indx = 3*(atom_index-1)
        if unit_lj:
            sigma = math.sqrt(temp/mass) # set mass unity; kB=1
        else: # real unit; kcal to kj; A/fs
            sigma=np.sqrt(4.184*temp/mass)/np.power(10,2)
        for j in range(3):
            v[indx+j] = random.gauss(0.0,sigma)
        mylmp.scatter_atoms("v",1,3,v)    # send it to LAMMPS
        mylmp.command("run 0 post no # temp=%f unit_lj?=%s velocity generated"%(temp,unit_lj))            # internal update

    ###################
    # trial insert
    ###################
    def insert_molecule(self,mylmp=None,dic_list=None,dic_selected=None,where=None,output=None):
        # dictionary
        dic_mc=dic_list['mc']
        dic_inputs=dic_list['inputs']
        #
        natoms_old=mylmp.get_natoms()
        if dic_mc['select_water']:  # water
            nstart=[natoms_old+x+1 for x in np.arange(dic_selected["natoms"])]
            dic_mc['selected']['water']=natoms_old+1    # oxygen index = the first atom
            m="insert water"
        elif dic_mc['select_lj']:  # lj
            nstart=[natoms_old+x+1 for x in np.arange(dic_selected["natoms"])]
            dic_mc['selected']['lj']=natoms_old+1    # oxygen index = the first atom
            m="insert lj"
        elif dic_mc['select_salt']: # salt
            nstart=[natoms_old+x+1 for x in np.arange(dic_selected["natoms"])] 
            dic_mc['selected']['salt']['cation']=natoms_old+1
            dic_mc['selected']['salt']['anion']=natoms_old+2
            m="insert salt-pair"
        # box info
        boxlo,boxhi,xy,yz,xz,periodicity,box_change=mylmp.extract_box()
        # random number for generation
        rseed1=random.randrange(1,100000)
        rseed2=random.randrange(1,100000)
        # water or lj
        if not dic_mc['select_salt']:
            if dic_selected["read_input"]:  # water
                # read molecule input file to insert single molecule
                mylmp.command("variable      offnatom    equal   %d   # offset" % 0)
                mylmp.command("create_atoms  ${offnatom} random  %d  %d  %s  mol %s %d units box # read" % (1, rseed1, dic_mc['region_to_insert'], dic_selected['mol_name'], rseed2))
                need_vel_gen=True
            else:   # lj solvent
                mylmp.command("create_atoms  %s random  %d  %d  %s units box # not_read" % (dic_selected['type_neq_s'], 1, rseed1, dic_mc['region_to_insert']))
                # velocity
                self.init_vel_single_lj(mylmp=mylmp,dic=dic_inputs['neqMD'],where="trial_insert",temp=dic_inputs['thermal'],atom_index=nstart[0],mass=dic_selected['mass'][0])
                ### test velocity
                v=mylmp.gather_atoms("v",1,3)
                for e,iontype in zip(nstart,["lj_solvent"]):
                    indx = 3*(e-1)
                    mylmp.command("# generated_velocity_lj %d %s vel %f %f %f"%(e,iontype,v[indx],v[indx+1],v[indx+2]))

                need_vel_gen=False
            # define a group of the created molecule
            mylmp.command("group   widom_test  type    %s # %s" % (dic_selected['type_neq_s'],m))
            # post processing
            if need_vel_gen:    # create a flyingn water
                # random orientation/position of water molecule
                oxy,hyd,imageh,vwat=mywater.generate_water(mylmp=mylmp,dic_inputs=dic_inputs,mass=dic_selected["mass"],boxlo=boxlo,boxhi=boxhi)
                # oxygen atom
                mylmp.command("set atom %d x %f y %f z %f"%(nstart[0],oxy[0],oxy[1],oxy[2]))
                mylmp.command("set atom %d image 0 0 0"%nstart[0])
                # hyrogen positions/images
                for j in np.arange(2):
                    mylmp.command("set atom %d x %f y %f z %f"%(nstart[1]+j,hyd[j][0],hyd[j][1],hyd[j][2]))
                    mylmp.command("set atom %d image %d %d %d"%(nstart[1]+j,imageh[j][0],imageh[j][1],imageh[j][2]))
                # velocities
                mylmp.command("set atom %d vx %f vy %f vz %f # water-O velocity"%(nstart[0],vwat[0],vwat[1],vwat[2]))
                mylmp.command("set atom %d vx %f vy %f vz %f # water-H velocity"%(nstart[1],vwat[3],vwat[4],vwat[5]))
                mylmp.command("set atom %d vx %f vy %f vz %f # water-H velocity"%(nstart[1]+1,vwat[6],vwat[7],vwat[8]))
            f_instant=False # no instant MC for water/lj
        else:   # salt-pair generation with bias
            # bias in interionic distance
            bias_f=dic_mc['bias']['bias_f']
            alpha_bias=dic_mc['bias']['alpha']
            x_center=dic_mc['bias']['center']
            rcut=dic_mc['bias']['ex_vol']
            mylmp.command("#The center (%f) and alpha (%f) of %s biasing function"%(x_center,alpha_bias,bias_f))
            # generate a flying ion-pair
            pos_c,pos_a,prob,image_a,f_instant=mysalt.generate_salt(mylmp=mylmp,boxlo=boxlo,boxhi=boxhi,bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center,rcut=rcut)
            dic_mc["bias"]["factor_ins"]=1./prob  # Volume if random insertion
            # positions of ions
            for e,ion in zip(dic_selected['type_neq_l'],[pos_c,pos_a]):
                mylmp.command("create_atoms %d single %f %f %f"%(e,ion[0],ion[1],ion[2]))
            # image for anion; due to pbc
            mylmp.command("set atom %d image %d %d %d"%(nstart[1],image_a[0],image_a[1],image_a[2]))
            # define group
            mylmp.command("group   widom_test  type    %s # %s" % (dic_selected['type_neq_s'],m))
            ## assign new momenta for flying ions
            #self.init_vel_single(mylmp=mylmp,dic=dic_inputs['neqMD'],where="trial_insert",temp=dic_inputs['temp'])
            for j,ion in enumerate(nstart):
                self.init_vel_single_lj(mylmp=mylmp,dic=dic_inputs['neqMD'],where="trial_insert",temp=dic_inputs['thermal'], atom_index=ion, mass=dic_selected['mass'][j], unit_lj=dic_mc['lj_flag'])
            #### test velocity
            v=mylmp.gather_atoms("v",1,3)
            for e,ion,iontype in zip(nstart,[pos_c,pos_a],["cation","anion"]):
                indx = 3*(e-1)
                mylmp.command("# generated_velocity_salt %d %s vel %f %f %f"%(e,iontype,v[indx],v[indx+1],v[indx+2]))
        ###
        # assign partial charges
        for t,c in zip(dic_selected['type_neq_l'],dic_selected['charge']):
            mylmp.command("set type %d charge %f" % (t,c))
        # internally save positions
        mylmp.command("run 0 post no # insert flying molecule")      # submit a new molecule internally
        natoms_new=mylmp.get_natoms()
        return natoms_new,"New flying particles assigned: (%d <- %d/%s)"%(natoms_new,natoms_old,where),f_instant

    def delete_molecule(self, mylmp=None, ljcut=9., where=None,ff='SD_SPCE'):
        if where == "accept_delete_lj":
            if ff != "simple_lj":
                mylmp.command("delete_atoms group widom_test compress no bond yes mol yes # accept_delete_lj")
                mylmp.command("run 0 post no # accept deleting flying water/salt molecule")    # submit a new particle internally
            else:
                mylmp.command("delete_atoms group widom_test # lj compress no # lj bond yes mol yes")
                mylmp.command("run 0 post no # accept deleting flying lj molecule")    # submit a new particle internally
            return "Delete flying particles: (%s)"%(where)
        if where != "reject_insert":
            # new method to make it ideal particles
            mylmp.command("# dealing with ideal flying particles: force field=%s"%ff)
            mylmp.commands_string(mypair.pair_coeff_delete_trial(ljcut=ljcut, ff=ff))
            mylmp.command("run 0 post no # making flying moelcule non-interacting")    # submit a new particle internally
            return "Delete flying particles (ideal gas): (%s)"%(where)
        else:
            # Remove the created molecule
            if ff != "simple_lj":
                mylmp.command("delete_atoms group widom_test compress no bond yes mol yes # %s"%where)
            else:
                mylmp.command("delete_atoms group widom_test compress no # lj bond yes mol yes # %s"%where)
                #mylmp.command("reset_atom_ids sort yes # lj")
            mylmp.command("run 0 post no # delete trial-inserted flying molecule")    # submit a new particle internally
            return "Delete flying particles: (%s)"%(where)

    def select_molecule(self,mylmp=None,dic_mc=None,dic_selected=None,x=None,where=None):
        # select delete molecules
        if dic_mc['select_lj']:         # for lj solvent
            ndelete=random.choice(dic_mc['prop_to_save']['lj_list'])
            selected=" {} ".format(ndelete)
            dic_mc['selected']['lj']=ndelete
            message="lj solvent is selected: %d\n"%(ndelete)
            #### test velocity
            v=mylmp.gather_atoms("v",1,3)
            for e,iontype in zip([ndelete],["lj-solvent"]):
                indx = 3*(e-1)
                mylmp.command("# selected_velocity_lj %d %s vel %f %f %f"%(e,iontype,v[indx],v[indx+1],v[indx+2]))
        elif dic_mc['select_water']:    # for water
            ndelete=random.randint(1,dic_mc['prop_to_save']['nwater'])  # random selection with oxygen atom index
            mylmp.command("# water selected (%d out of %d)"%(ndelete,len(dic_mc['prop_to_save']['water_atom_list'])))
            water_atom_list=dic_mc['prop_to_save']['water_atom_list'][ndelete-1]
            water_atom_list.sort(key=lambda x: x[1]) # from oxygen to hydrogens
            # based on the data structure - not in use; since the data is not consecutive in order
            #ndelete=random.choice(dic_mc['prop_to_save']['water_list'])
            #selected=" {} {} {}".format(ndelete,ndelete+1,ndelete+2)
            selected=" {} {} {}".format(water_atom_list[0][0],water_atom_list[1][0],water_atom_list[2][0])
            dic_mc['selected']['water']=ndelete
            message="water is selected: %d %d %d\n"%(water_atom_list[0][0],water_atom_list[1][0],water_atom_list[2][0])
        elif dic_mc['select_salt']:     # for salt
            # bias in interionic distance
            bias_f=dic_mc['bias']['bias_f']
            alpha_bias=dic_mc['bias']['alpha']
            x_center=dic_mc['bias']['center']
            rcut=dic_mc['bias']['ex_vol']
            # random choice of a cation and its 3D position
            ndelete_c=random.choice(dic_mc['prop_to_save']['cations_list'])
            atom_c=3*(ndelete_c-1)
            pos_c=[x[atom_c],x[atom_c+1],x[atom_c+2]]
            mylmp.command("#chosen_cation deletion_move %f %f %f # select_molecule"%(x[atom_c],x[atom_c+1],x[atom_c+2]))
            # bias choice of anion
            pw=0.
            choice_anion_prob=[]
            for ndelete_a in dic_mc['prop_to_save']['anions_list']:
                atom_a=3*(ndelete_a-1)
                tri_a=[x[atom_a],x[atom_a+1],x[atom_a+2]]
                p=mysalt.bias_salt_del_choice(mylmp=mylmp,salt=[pos_c,tri_a],bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                pw+=p
                choice_anion_prob.append(pw)
            # chosen anion
            choice_anion_prob=[x/pw for x in choice_anion_prob]
            rand=random.random()
            for p,ndelete_a in zip(choice_anion_prob,dic_mc['prop_to_save']['anions_list']):
                atom_a=3*(ndelete_a-1)
                pos_a=[x[atom_a],x[atom_a+1],x[atom_a+2]]
                if rand<p:
                    prob=mysalt.bias_salt_del_choice(mylmp=mylmp,salt=[pos_c,pos_a],bias_f=bias_f,alpha_bias=alpha_bias,x_center=x_center)
                    mylmp.command("#chosen_anion deletion_move %f %f %f # select_molecule"%(x[atom_a],x[atom_a+1],x[atom_a+2]))
                    break
            dic_mc['bias']["factor_del"]=pw/prob    # N+1 if equally probable
            selected=" {} {}".format(ndelete_c,ndelete_a)
            dic_mc['selected']['salt']['cation']=ndelete_c
            dic_mc['selected']['salt']['anion']=ndelete_a
            message="salt is selected: %d %d\n"%(ndelete_c,ndelete_a)
            #### test velocity
            v=mylmp.gather_atoms("v",1,3)
            for e,ion,iontype in zip([ndelete_c,ndelete_a],[pos_c,pos_a],["cation","anion"]):
                indx = 3*(e-1)
                mylmp.command("# selected_velocity_salt %d %s vel %f %f %f"%(e,iontype,v[indx],v[indx+1],v[indx+2]))
        # define group of widom_test
        mylmp.command("group widom_test id %s # select_molecule" % selected)
        # in case of trial delete move is rejected - match with full_list
        dic_mc['selected']['list']=[int(x) for x in selected.split()]
        # switch type from eq to neq
        self.type_change_molecule(mylmp=mylmp,dic_mc=dic_mc,dic_selected=dic_selected,where=where)
        mylmp.command("run 0 post no # done selecting a flying molecule")    # submit a new particle internally
        return message+"Select flying particles (%s): (%s)"%(selected,where)

    def update_list_atoms(self,mylmp=None,dic_mc=None,output=None, ordering=True):
        # initialize arrays
        dic_mc['prop_to_save']['lj_list']=[]
        dic_mc['prop_to_save']['water_list']=[]
        dic_mc['prop_to_save']['cations_list']=[]
        dic_mc['prop_to_save']['anions_list']=[]
        # gather type
        natoms=mylmp.get_natoms()
        atom_id_new=(3*natoms*c_int)()
        atom_type = mylmp.gather_atoms("type",0,1)
        atom_id = mylmp.gather_atoms("id",0,1)
        atom_molid = mylmp.gather_atoms("molecule",0,1)
        max_molid=max(atom_molid)
        mylmp.command("# updating atom list by mol_id %d"%max_molid)
        if max_molid<1:
            mylmp.command("# this is a pure atomic-system atom_id (%d)"%max_molid)
        else:
            # need to redo numbering according to molid
            mol_list=[[] for x in np.arange(max_molid)]
            for idx,at in enumerate(atom_type):
                mid=atom_molid[idx] # molecule id
                aid=atom_id[idx]    # atom id
                mol_list[mid-1].append([aid,at])
            for m,moll in enumerate(mol_list):
                if len(moll)<2: # atomic molecule
                    mol_list.remove(moll) 
                else:
                    moll.sort() # by default, take the first parameter
            dic_mc['prop_to_save']['water_atom_list']=mol_list
            mylmp.command("# Number of molecules (eg water) in updating index: %d"%len(mol_list))
            del mol_list
        # index ordering
        for idx,at in enumerate(atom_type):
            aid=atom_id[idx]    # atom id
            if at in dic_mc["atom_type"]["lj"]:
                dic_mc['prop_to_save']['lj_list'].append(aid)
            elif at in dic_mc["atom_type"]["water"]:
                dic_mc['prop_to_save']['water_list'].append(aid)
            elif at in dic_mc["atom_type"]["cation"]:
                dic_mc['prop_to_save']['cations_list'].append(aid)
            elif at in dic_mc["atom_type"]["anion"]:
                dic_mc['prop_to_save']['anions_list'].append(aid)
            #else: # for instance hydrogen in water
            #    print ("Error in updating atom list")
            #    exit()
        dic_mc['prop_to_save']['nlj']=len(dic_mc['prop_to_save']['lj_list'])
        dic_mc['prop_to_save']['nwater']=len(dic_mc['prop_to_save']['water_list'])
        dic_mc['prop_to_save']['nsalt']=len(dic_mc['prop_to_save']['cations_list'])
        print ("Number of lj-solvent/water/salt-pair: %d/%d/%d"%(dic_mc['prop_to_save']['nlj'],dic_mc['prop_to_save']['nwater'],dic_mc['prop_to_save']['nsalt']),file=output)
        if dic_mc['prop_to_save']['nsalt'] != len(dic_mc['prop_to_save']['anions_list']):
            print ("ERROR (update_list_atoms): not consistent number of salt")
            print ("cation/anion",len(dic_mc['prop_to_save']['cations_list']),len(dic_mc['prop_to_save']['anions_list']))
            exit()
        mylmp.command("# done updating numbering after accepted exchange")
        #if select_lj:
        #    #last=self.dic_mc['selected']['lj']['solvent']
        #    #self.dic_mc['prop_to_save']['lj_list_solvent'].remove(last)
        #    self.dic_mc['prop_to_save']['nlj_solvent']-=1
        #else:
        #    #last=self.dic_mc['selected']['lj']['solute']
        #    #self.dic_mc['prop_to_save']['lj_list_solute'].remove(last)
        #    self.dic_mc['prop_to_save']['nlj_solute']-=1


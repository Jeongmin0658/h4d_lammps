from ctypes import *
import sys,random,math
import numpy as np
import os

#me = 0
#from mpi4py import MPI
#me = MPI.COMM_WORLD.Get_rank()
#nprocs = MPI.COMM_WORLD.Get_size()

# import modules
from dictionaries import LMP_dics
mydics=LMP_dics()
from define_lmp_misc import LMP_misc
mymisc=LMP_misc()
from define_trial_moves import LMP_trial_moves
mymove=LMP_trial_moves()
from define_bias_salt import LMP_define_bias_salt
mysalt=LMP_define_bias_salt()
# selection process
from define_select_move import LMP_define_select_move
myselect=LMP_define_select_move()

# Class for hybrid MD/MC
class LMP_MC():
    def __init__(self,cut_sel=12,width_sel=0.1,selection=False):
        self.cut_sel=cut_sel 
        self.width_sel=width_sel
        self.selection=selection
        self.intro="Define Monte Carlo (Selection?%s/cut=%f/width=%f)"%(self.selection,self.cut_sel,self.width_sel)

    def __clear(self):
        del self.dic_selected
        self.dic_mc['bias']['selection']['ins']=1.
        self.dic_mc['bias']['selection']['del']=1.
        
    # main driver to run MC
    def run_MC(self, mylmp=None,output="mc_test",dic_list=None,flag_instant=False,flag_neq=False):
        self.flag_neq=flag_neq
        self.dic_inputs=dic_list['inputs']
        self.dic_mc=dic_list['mc']
        self.dic_flying=dic_list['flying']
        dic_peratom=dic_list['peratom']
        # update current status
        mydics.update_peratom(mylmp=mylmp,status="current",dic_peratom=self.dic_peratom)
        if self.selection:
            # selection bias factor after NEMD
            fact_sel=myselect.get_factor_selection(mylmp=mylmp, dic_mc=self.dic_mc, output=output, cut=self.cut_sel, width=self.width_sel)
            mylmp.command("# selection bias after (%f/%f) = %f"%(self.cut_sel,self.width_sel,fact_sel))
            if self.dic_mc['select_insert']:
                self.dic_mc['bias']['selection']['del']=fact_sel
            else:
                self.dic_mc['bias']['selection']['ins']=fact_sel
            # calculate bias factor after NEMD
        if self.dic_mc['select_salt']:    # salt
            mysalt.bias_after_NEMD(mylmp=mylmp,x=dic_peratom["current"]["x"],forward=self.dic_mc['select_insert'],dic_mc=self.dic_mc)
        #####################
        # Metropolis MC
        #####################
        accept,dke,dpe,dE,m=self.check_accept(mylmp=mylmp, output=output)

        # need to write after calculating bias 
        if self.dic_mc['select_salt']:
            print ("%s (dke=%f; dpe=%f; dE=%f; vol=%f; prob_ins_inv=%e; prob_del=%e); intant_MC=%s; not_to_sample=%s; @ RT=%f"%(m,dke,dpe,dE,self.dic_mc['volume'],self.dic_mc['bias']['factor_ins'],self.dic_mc['bias']['factor_del'],flag_instant,flag_neq,self.dic_inputs["thermal"]),file=output)
        elif self.dic_mc['select_lj']:
            print ("%s (dke=%f; dpe=%f; dE=%f; vol=%f; prob_ins_inv=%e; prob_del=%e); intant_MC=%s; not_to_sample=%s; @ RT=%f"%(m,dke,dpe,dE,self.dic_mc['volume'], self.dic_mc['volume'], self.dic_mc['prop_to_save']['nlj'], flag_instant,flag_neq,self.dic_inputs["thermal"]),file=output)
        else:
            print ("%s (dke=%f; dpe=%f; dE=%f; vol=%f; prob_ins_inv=%e; prob_del=%e); intant_MC=%s; not_to_sample=%s; @ RT=%f"%(m,dke,dpe,dE,self.dic_mc['volume'], self.dic_mc['volume'], self.dic_mc['prop_to_save']['nwater'], flag_instant,flag_neq,self.dic_inputs["thermal"]),file=output)
        #####################
        # post processing after MC acceptance
        #####################
        m=self.post_check_accept(mylmp=mylmp,dic_selected=self.dic_selected,output=output)
        print (m,file=output)
        # clear step
        self.__clear()

    def save_old_state(self,mylmp=None, output="mc_test",dic_list=None):
        # dictionary
        self.dic_inputs=dic_list['inputs']
        self.dic_mc=dic_list['mc']
        self.dic_peratom=dic_list['peratom']
        self.dic_flying=dic_list['flying']

        # initialize peratom array
        mydics.init_peratom(mylmp=mylmp,dic_peratom=self.dic_peratom)
        # save OLD state and pe_old/ke_old
        natoms_old=mylmp.extract_global("natoms",0)
        # update peratom array
        mydics.update_peratom(mylmp=mylmp,status="old",dic_peratom=self.dic_peratom)
        # save OLD state and pe_old/ke_old
        natoms_old=mylmp.extract_global("natoms",0)
        # note: kinetic energy will be calculated after creating/selecting flying particles
        pe_old=mylmp.extract_compute("thermo_pe",0,0)
        # compute intramolecular energies
        bond_old,angle_old=mymisc.compute_water_intramolecular(mylmp=mylmp)
        self.dic_inputs['neqMD']['energies']['bond_w_old']=bond_old
        self.dic_inputs['neqMD']['energies']['angle_w_old']=angle_old
        mol_old=bond_old+angle_old
        pe_old-=mol_old # no need since eqMD evolves with shake or rigid/nve
        self.dic_inputs['neqMD']['energies']['pe_old']=pe_old
        ke_old=mylmp.extract_compute("thermo_ke",0,0) # dummy to check
        vol_old=mymisc.get_volume(mylmp=mylmp) # volume
        print ("old state right after eqMD (natom=%d,pe_old=%f,ke_old=%f,bond+angle=%f,vol=%f)"%(natoms_old,pe_old,ke_old,mol_old,vol_old),file=output)

    def prep_MC(self,mylmp=None, output="mc_test",dic_list=None,istep=0):
        # dictionary
        self.dic_inputs=dic_list['inputs']
        self.dic_mc=dic_list['mc']
        self.dic_peratom=dic_list['peratom']
        self.dic_flying=dic_list['flying']

        ## initialize peratom array
        #mydics.init_peratom(mylmp=mylmp,dic_peratom=self.dic_peratom)
        ## update peratom array
        #mydics.update_peratom(mylmp=mylmp,status="old",dic_peratom=self.dic_peratom)
        # save OLD state and pe_old/ke_old
        natoms_old=mylmp.extract_global("natoms",0)
        ## note: kinetic energy will be calculated after creating/selecting flying particles
        #pe_old=mylmp.extract_compute("thermo_pe",0,0)
        ## compute intramolecular energies
        #bond_old,angle_old=mymisc.compute_water_intramolecular(mylmp=mylmp)
        #self.dic_inputs['neqMD']['energies']['bond_w_old']=bond_old
        #self.dic_inputs['neqMD']['energies']['angle_w_old']=angle_old
        #mol_old=bond_old+angle_old
        #pe_old-=mol_old # no need since eqMD evolves with shake or rigid/nve
        #self.dic_inputs['neqMD']['energies']['pe_old']=pe_old
        #ke_old=mylmp.extract_compute("thermo_ke",0,0) # dummy to check
        #vol_old=mymisc.get_volume(mylmp=mylmp) # volume
        #print ("old state right after eqMD (natom=%d,pe_old=%f,ke_old=%f,bond+angle=%f,vol=%f)"%(natoms_old,pe_old,ke_old,mol_old,vol_old),file=output)

        # probabilities
        prob_wat=self.dic_mc['prob_water']
        prob_lj=self.dic_mc['prob_lj']+prob_wat
        prob_salt=self.dic_mc['prob_salt']+prob_lj
        # insert or delete
        prob_insert=self.dic_mc['prob_insert']
        ###
        # trial move
        ###
        self.dic_mc['select_water']=False
        self.dic_mc['select_lj']=False
        self.dic_mc['select_salt']=False
        # randomly select a flying molecule
        rand=random.random()
        if rand < prob_wat:     # water
            self.dic_mc['select_water']=True
            self.dic_selected=self.dic_flying['water']
            print ("water is picked", self.dic_mc['select_water'],file=output)
        elif rand < prob_lj: # lj
            self.dic_mc['select_lj']=True
            self.dic_selected=self.dic_flying['lj']
            print ("lj is picked", self.dic_mc['select_lj'],file=output)
        else:                   # salt
            self.dic_mc['select_salt']=True
            self.dic_selected=self.dic_flying['salt']
            print ("salt is picked", self.dic_mc['select_salt'],file=output)
        # pick a move
        rand=random.random()
        # check if this would be an instant MC
        f_instant=False # false = hybrid NEMD/MC; you just choose
        if rand<prob_insert:
            self.dic_mc['select_insert']=True
            self.dic_inputs['neqMD']['neq']['insertion']=True
            # create flying particles to insert and assign new positions
            natoms_new,m,f_instant=mymove.insert_molecule(mylmp=mylmp,dic_list=dic_list,dic_selected=self.dic_selected,where="trial_insert",output=output)
            # update peratom dic
            mydics.update_peratom(mylmp=mylmp,status="current",dic_peratom=self.dic_peratom)
            print ("%s (%d <= %d)"%(m,natoms_new,natoms_old),file=output)
        else:
            self.dic_mc['select_insert']=False
            self.dic_inputs['neqMD']['neq']['insertion']=False
            # select particles to delete with dic_selected dictionary
            m=mymove.select_molecule(mylmp=mylmp,dic_mc=self.dic_mc,dic_selected=self.dic_selected,x=self.dic_peratom["old"]["x"],where="trial_delete")
            print (m,file=output)
        ###
        # evolving groups - restricted region around flying ions
        if self.selection:
            myselect.select_free_to_move(mylmp=mylmp, dic_mc=self.dic_mc, output=output, cut=self.cut_sel, width=self.width_sel)
            fact_sel=myselect.get_factor_selection(mylmp=mylmp, dic_mc=self.dic_mc, output=output, cut=self.cut_sel, width=self.width_sel)
            mylmp.command("# selection bias before (%f/%f) = %f"%(self.cut_sel,self.width_sel,fact_sel))
            if self.dic_mc['select_insert']:
                self.dic_mc['bias']['selection']['ins']=fact_sel
            else:
                self.dic_mc['bias']['selection']['del']=fact_sel
            # group setting for non-rigid lj or salt
            if self.dic_mc['select_salt'] or self.dic_mc['select_lj']:
                mylmp.command("group noneq_run          subtract     noneq_free_other widom_test # %s"%self.dic_inputs['neqMD']['evolve']['flying'])
                mylmp.command("group noneq_run_rigid    intersect   water            noneq_free_water")
            # group setting for rigid water
            # union need to be changed to be subtract when the selection is in use
            elif self.dic_mc['select_water']:
                mylmp.command("group noneq_run          subtract     noneq_free_other widom_test # %s"%self.dic_inputs['neqMD']['evolve']['flying'])
                mylmp.command("group water_to_relax     intersect   water            noneq_free_water")
                mylmp.command("# comment: union should be subtract when the selection process is in use")
                mylmp.command("group noneq_run_rigid    union       water_to_relax   widom_test # %s"%self.dic_inputs['neqMD']['evolve']['flying'])
                mylmp.command("group water_to_relax     delete")
            mylmp.command("#group noneq_run_rigid        subtract    noneq_run_rigid_pre widom_test")
            mylmp.command("#group noneq_run_rigid_pre    delete")
            mylmp.command("group noneq_free_other       delete")
            mylmp.command("group noneq_free_water       delete")
            mylmp.command("run 0 post no # after selection of all moving atoms")
        else:
           # group setting for non-rigid lj or salt
            if self.dic_mc['select_salt']:
                mylmp.command("group noneq_run_rigid    intersect   all         water               # select_salt no selection")
                mylmp.command("group noneq_run_pre      subtract    all         noneq_run_rigid     # select_salt no selection")
                mylmp.command("group noneq_run          subtract    noneq_run_pre   widom_test      # select_salt no selection")
                mylmp.command("group noneq_run_pre      delete")
            elif self.dic_mc['select_lj']:
                mylmp.command("group noneq_run_rigid    intersect   all         water               # select_lj no selection")
                mylmp.command("group noneq_run          subtract    all         noneq_run_rigid     # select_lj no selection")
            # group setting for rigid water
            # union need to be changed to be subtract when the selection is in use
            elif self.dic_mc['select_water']:
                mylmp.command("group noneq_run_rigid    union       water       widom_test          # select_water no selection")
                mylmp.command("group noneq_run          subtract    all         noneq_run_rigid     # select_water no selection") 
        ###
        # total momentum after assigning new velocities for flying particles
        ke_old = mylmp.extract_compute("thermo_ke",0,0)
        pe_old = mylmp.extract_compute("thermo_pe",0,0) # dummy to check
        self.dic_inputs['neqMD']['energies']['ke_old']=ke_old
        print ("old state right before neqMD (natom=%d,pe_old=%f,ke_old=%f)"%(natoms_old,pe_old,ke_old),file=output)
        return f_instant

    # check if accept a trial move
    def check_accept(self,mylmp=None,output='mc_test'):
        # energies
        pe_new=self.dic_inputs['neqMD']['energies']['pe_new']
        ke_new=self.dic_inputs['neqMD']['energies']['ke_new']
        pe_old=self.dic_inputs['neqMD']['energies']['pe_old']
        ke_old=self.dic_inputs['neqMD']['energies']['ke_old']
        print ("Energy check:",pe_new,pe_old,ke_new,ke_old,file=output)
        # compute total energy differences during neMD for insertion/deletion
        dke=ke_new-ke_old           # kinetic energy difference
        dpe=pe_new-pe_old           # potential energy difference
        dE=dke+dpe                  # total energy difference
        # test acceptance
        accept,m=self.accept_rule(dE=dE,kT=self.dic_inputs['thermal'], output=output)
        return accept,dke,dpe,dE,"Acceptance tested: "+m

    def post_check_accept(self, mylmp=None, dic_selected=None, output="mc_test"):
        accept=self.dic_mc['accept']
        forward=self.dic_mc['select_insert']
        if not accept:  # reject = return to the previous state
            if forward:
                # delete the inserted flying molecule
                mymove.delete_molecule(mylmp=mylmp,ljcut=self.dic_mc['cutoff']['lj'],where="reject_insert",ff=self.dic_inputs['neqMD']['ff'])
                m="reject_insert:::delete the inserted flying molecule\n"
            else:
                # return the types
                mymove.accept_insert(mylmp=mylmp,dic_selected=dic_selected,where="reject_delete")
                # return the charges
                mymove.charge_return(mylmp=mylmp,dic_mc=self.dic_mc,dic_selected=dic_selected,where="reject_delete")
                m="reject_delete:::return the types so make them interacting again\n"
            # return to the old state
            mymove.return_position(mylmp=mylmp,dic_peratom=self.dic_peratom)
            m+="return to the old state using dic_peratom\n"
        else: # accepted move
            if forward:
                # change the type
                mymove.accept_insert(mylmp=mylmp,dic_selected=dic_selected,where="accept_insert")
                ## change the charge ?
                #mymove.charge_return(mylmp=mylmp,dic_mc=self.dic_mc,dic_selected=dic_selected,where="accept_insert")
                m="accept_insert:::change the type so they are no long flying ones and update lists\n"
            else:
                # delete the group
                mymove.delete_molecule(mylmp=mylmp, where="accept_delete_lj",ff=self.dic_inputs['neqMD']['ff'])
                m="accept_delete:::flying molecule is deleted\n"
            # contiguous numbering after the change in particle number
            mylmp.command("reset_atom_ids sort yes # id numbering")
            if dic_selected['molecule']:
                mylmp.command("reset_mol_ids all compress yes #single yes # molecule")
            mylmp.command("run 0 post no # just for ids re-setting")
            # update list of solvent and solute
            mymove.update_list_atoms(mylmp=mylmp,dic_mc=self.dic_mc,output=output, ordering=True)
        # delete widom_test group
        mylmp.command("group widom_test         delete")
        mylmp.command("group noneq_run          delete")
        mylmp.command("group noneq_run_rigid    delete")
        if self.selection:
            mylmp.command("group noneq_free         delete # after calculating selection factor")
        m+="delete widom_test group"
        return m

    # Metropolis acceptance criterion
    def accept_rule(self,dE=0.,kT=1., output=None):
        # reduced energy
        dE_red=dE/kT
        # compute the Boltzmann factor
        try:
            boltz=math.exp(-dE_red)
            m=""
        except OverflowError:
            if dE<0.:
                boltz = float('inf')
                m="WARNING:::TOO LARGE BOLTZ %f (dE=%f @ kT=%f)"%(boltz,dE,kT)
            else:
                return False,"ERROR in boltz (%f) with dE (%f)"%(boltz,dE)
        # gcmc
        if self.dic_mc['select_lj']:
            numb=self.dic_mc['prop_to_save']['nlj']     # number is not yet updated so +1 below
            activity=self.dic_mc['activity']['lj']
            f_bias=1.
        elif self.dic_mc['select_water']:
            numb=self.dic_mc['prop_to_save']['nwater']  # number is not yet updated so +1 below
            activity=self.dic_mc['activity']['water']
            f_bias=1.
        elif self.dic_mc['select_salt']:
            numb=self.dic_mc['prop_to_save']['nsalt']   
            activity=self.dic_mc['activity']['salt']
            # factor for the bias
            numb_bias=self.dic_mc['bias']['factor_del']   
            vol_bias=self.dic_mc['bias']['factor_ins']
            f_bias=vol_bias/numb_bias
        else:
            print ("ERROR:::wrong select_flag in accept_rule")
            exit()
        # calculate the factor for Metropolis criterion
        if self.dic_mc['select_insert']:     # insertion
            self.dic_mc['trial']['insert']+=1
            bias_sel=self.dic_mc['bias']['selection']['del']/self.dic_mc['bias']['selection']['ins']    # selection bias factor
            prep=self.dic_mc['volume']/float(numb+1)*f_bias*bias_sel                                    # include both bias factors
            try:
                factor=prep*math.exp(-dE_red+activity)  # for MC
                m="therm="+str(kT)+" (insertion)"
            except OverflowError:
                if dE_red-activity<0.:
                    factor = 1.0
                    m="WARNING:::TOO LARGE BOLTZ (insertion factor = 1) %f (dE=%f @ kT=%f)"%(dE,kT)
                else:
                    return False,"ERROR in boltz (%f) with dE (%f)"%(boltz,dE)
        else:                               # deletion - MC with non-interacting a flying molecule
            self.dic_mc['trial']['delete']+=1
            bias_sel=self.dic_mc['bias']['selection']['ins']/self.dic_mc['bias']['selection']['del']    # selection bias factor
            prep=float(numb)/self.dic_mc['volume']/f_bias*bias_sel                                      # include both bias factors
            try:
                factor=prep*math.exp(-dE_red-activity)  # for MC
                m="therm="+str(kT)+" (deletion)"
            except OverflowError:
                if dE_red+activity<0.:
                    factor = 1.0
                    m="WARNING:::TOO LARGE BOLTZ (deletion factor = 1) %f (dE=%f @ kT=%f)"%(dE,kT)
                else:
                    return False,"ERROR in boltz (%f) with dE (%f)"%(boltz,dE)
        print ("# bias selection to relax (accept_rule): %f\n# bias interionic distance(accept_rule): %f\n# prep = %e"%(bias_sel,f_bias,prep), file=output)
        # check if this is gcmc
        if not self.dic_mc['gcmc'] or self.dic_mc['loop']:
            self.dic_mc['accept']=False
            return False, "not_gcmc"
        if self.flag_neq:   # failed trial move; error in MD so not accept anyway
            return False, "gcmc_not_to_sample"
        ## TEST
        #self.dic_mc['accept']=True
        #return True, "test all accepted"
        if random.random() <= factor:
            self.dic_mc['n_accept']+=1
            self.dic_mc['accept']=True
            return True, m
        self.dic_mc['accept']=False
        return False, m 


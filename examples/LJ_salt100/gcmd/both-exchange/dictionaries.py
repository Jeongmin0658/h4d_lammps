import numpy as np
import math
import os,sys

#############################################################
### CHANGE HERE ###
# where your source folder is
source_where="WHERE_YOUR_SOURCES_ARE"
sys.path.append(source_where)
#############################################################

from define_lmp_misc import LMP_misc
mymisc=LMP_misc()

class LMP_dics():
    def __init__(self):
        """
        GCMD for bulk electrolytes
        Note: 
            + "None" will be updated during simulations
        """
        #############################################################
        # Preliminary setting
        #############################################################
        kJ_to_kcal=1.0/4.184
        ###
        # chemical potential for GCMD
        # in unit of kT 
        # De Broglie wavelength = 1
        chem_pot_lj=-1.754          # monoatomic LJ solvent
        chem_pot_water=-15.37       # SPC/E water (rigid)
        chem_pot_salt=-14.81        # a pair of monoatomic cation and anion
        ###

        # region to insert
        region_to_insert="NULL"     # Do NOT change this

        # force-field (See define_pair.py)
        ff='simple_lj'              # example:'SD_SPCE' # or 'JD_SPCE' or 'simple_lj'

        # unit in lammps
        lj_unit=True                # true if lammps uses lj unit, otherwise False

        # timestep for time evolution
        dt = 0.005                  # for equilibrium MD
        dt_neq = 0.01               # for NEMD

        # temperature/pressure
        temp=1.0 
        press=1.0

        # damping for equilibrium MD
        # note: NEMD will be in NVE ensemble
        temp_damp=100.0*dt          # for thermostat
        press_damp=1000.0*dt        # for barostat

        # LJ parameter setting
        lj_mixing='arithmetic'
        lj_tail='no'
        lj_shift='yes'

        # cutoff distances for LJ/Coulomb interactions
        cutoff_lj = 2.5
        cutoff_coul = 3.5
        cutoff_coul_neq = 3.5

        # implicit confinement
        # Do NOT change this; this is NOT used
        confinement = {
            'flag': False, #True,
            'wall': "wall/lj93",
            'eps': 8.0,
            'sig': 1.0,
            'cut': 40.0,
            'spacing': 1.0/math.sqrt(2.0)
        }

        #############################################################
        # MC inputs
        #############################################################
        self.dic_mc={
            "ref"           : "Luc Belloni, J. Chem. Phys. 149, 094111 (2018)",
            "gcmc"          : True,         # if False, no trial move is accepted, but measure chemical potential
            "loop"          : False,        # NOT in use
            "lj_flag"       : lj_unit,      # lj units
            "atom_type"     : {             # Consistent with your lammps input
                "lj"        : [1],
                "water"     : [],          
                "cation"    : [2],
                "anion"     : [3]
            },
            "bias"          :{              # for ion-pair exchange
                "bias_f"    : "bimodal",    # either bimodal, gaussian, uniform
                "alpha"     : 5.0,          # inverse of width (alpha_b in the paper)
                "center"    : 0.8,          # two center of bimodal (x_b in the paper)
                "ex_vol"    : 0.3,          # excluded volume (length) (ex_vol^3 = V_ex in the paper)
                "factor_ins"  : None,       # factor in acceptance prob. with the bias (B^del/B^ins)
                "factor_del"  : None,       # factor in acceptance prob. with the bias (B^ins/B^del)
                "selection": {              # NOT used in this version
                    "ins": 1.,
                    "del": 1.
                }
            },
            # probability of selecting species for flying particles
            # water + lj + salt = 1
            # go down "evolve" in neqMD section, and change it accordingly
            # example: "substract" for salt-exchange; "union" for water/lj exchange
            "prob_water"    : 0.0,           # electrolytes
            "prob_lj"       : 0.8,          # simple lj solvent
            "prob_salt"     : 0.2,          # electrolytes
            # probability of insertion (p_del = 1 - p_ins)
            "prob_insert"   : 0.5,          # for all ff; either 0.5 or 1 or 0
            # select: will be updated
            "select_insert" : None,
            "select_lj"     : None,         
            "select_water"  : None,
            "select_salt"   : None,
            "selected"      : {
                "list"      : None,
                "lj"        : None,
                "water"     : None,
                "salt"      :{
                    "cation": None,
                    "anion" : None
                },
            },
            "trial"         : {
                "insert": 0,
                "delete": 0
            },
            "confined"      : confinement.copy()["flag"],   # NOT used in this version
            "region_to_insert": region_to_insert,
            "accept"        : None,
            "n_accept"      : 0,
            "activity"      : {
                "lj"        : None,         # simple_lj
                "water"     : None,
                "salt"      : None
            },
            "volume"        : None,
            "cutoff"        : {             # will be updated
                "lj"        : None, 
                "coul"      : None
            },
            "prop_to_save"  : { 
                "ex_chem"   : {
                    "insert": 0.,
                    "delete": 0.,
                    "list"  :{
                        "insert": [],
                        "delete": []
                    }
                },
                "nwater"        : 0,
                "nsalt"         : 0,
                "nlj"           : 0,
                "lj_list"       : None,
                "water_list"    : None,
                "water_atom_list"    : None,
                "cations_list"  : None,
                "anions_list"   : None
            }
        }

        #############################################################
        # Flying molecule inputs
        # Be consistent with your lammps input
        #############################################################
        self.dic_flying={
            "lj":{ 
                "molecule"  : False, #True,
                "read_input": False, # if True, need to specify how to generate velocity
                "region_to_insert": region_to_insert,
                "natoms"    : 1,
                "group"     : ["ljs"],
                "type_eq_s" : "1",
                "type_eq_l" : [1],
                "type_neq_s": "6",      # string
                "type_neq_l": [6],      # list
                "mass"      : [1.0],  
                "charge"    : [0.],
                "full_list" : { 
                    "charge": [0.],
                    "type_eq_l" : [1],
                    "type_neq_l": [6],
                },
                "file"      : None, 
                "mol_name"  : "widom_lj",
                "activity"  : None      # will be updated according to chemical potential
            },
            "salt":{ 
                "molecule"  : False, 
                "read_input": False,    # if True, need to specify how to generate velocity
                "region_to_insert": region_to_insert,
                "natoms"    : 2,
                "group"     : ["ljs"],
                "type_eq_s" : "2 3",
                "type_eq_l" : [2,3],
                "type_neq_s": "7 8",    # string
                "type_neq_l": [7,8],    # list
                "mass"      : [1.0, 1.0],
                "charge"    : [1.0,-1.0],
                "full_list" : { # full list (match with dic_mc[selected][list])
                    "charge": [1.0,-1.0],
                    "type_neq_l": [7,8],
                    "type_eq_l" : [2,3]
                },
                "file"      : None, 
                "mol_name"  : "widom_lj",
                "activity"  : None 
            },
            "water":{
                "molecule"  : True, # SPC/E water consists of more than one atoms
                "read_input": True, # see "file" below
                "region_to_insert": region_to_insert,
                "natoms"    : 3,
                "group"     : ["oxygens","hydrogens"],
                "type_eq_s" : "1 2",
                "type_eq_l" : [1, 2],
                "type_neq_s": "8 9",    # string
                "type_neq_l": [8, 9],   # list
                "mass"      : [15.999400, 1.007940, 1.007940],
                "charge"    : [-0.847600,0.423800],
                "full_list" : { # full list (match with dic_mc[selected][list])
                    "charge": [-0.847600,0.423800,0.423800],
                    "type_neq_l": [8,9,9],
                    "type_eq_l" : [1,2,2]
                },
                "file"      : source_where+"/input_lmp/water_spce_widom.txt",
                "mol_name"  : "widom_water",
                "activity"  : None 
            }
        }

        # update chemical potential
        self.dic_flying['water']['activity']=chem_pot_water
        self.dic_flying['salt']['activity']=chem_pot_salt
        self.dic_flying['lj']['activity']=chem_pot_lj
        
        #############################################################
        # MD input files
        #############################################################
        self.dic_inputs={
            'random_seed':  9071,
            'input_file':{  # input lammps file to read
                'to_read':  'lmp_spm_bulk_re',
                'where':    './'
            },
            'output_file':  'output_test',  # NOT used in this version
            'confinement_factor': 1.0,      # NOT used in this version
            'cutoff':{
                'lj':       cutoff_lj,
                'coul':     cutoff_coul
            },
            'intramol':{    # for SPC/E water
                'bond':     'harmonic',
                'angle':    'harmonic'
            },
            'LJ_mixing':    lj_mixing,
            'LJ_tail':      lj_tail,
            'LJ_shift':     lj_shift,
            'thermal': None,            # will be updated
            'temp':  temp,
            'temp_damp': temp_damp,       # langevin damping parameter
            #########################
            # pre MD before GCMD
            #########################
            'eqMD-pre':{
                'mol': False,
                'ff': ff,
                'neq':{
                    'flag': False
                },
                'dt': dt,
                'temp': temp,
                'temp_damp': temp_damp,   # langevin
                'press': press,
                'press_damp': press_damp,
                'neigh_freq': 10,
                'thermo_freq': 4000,
                'style': 'lj/cut/coul/long',
                'ensemble': 'nvt',
                'dump': {
                    'flag': False,
                    'atom': "all",
                    'freq': 1000,
                    'file': "spce_bulk_neq.lammpstrj",
                    'style': "id type xu yu zu #vx vy vz"
                },
                'confinement': confinement.copy(),
                'fixid': ['eqmd','eql'],
                'cutoff':{
                    'lj':   cutoff_lj,
                    'coul': cutoff_coul
                },
                'intramol':{
                    'bond': 'harmonic',
                    'angle': 'harmonic'
                },
                'LJ_mixing' : lj_mixing,
                'LJ_tail'   : lj_tail,
                'LJ_shift'  : lj_shift,
                'shake': {
                    'flag'  :   False,
                    'fix'   :   "fix fOHs        all    shake 0.0001 20 ${rattle_freq} b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
            #########################
            # Equilibrium MD
            #########################
            'eqMD':{
                'mol': False,
                'ff': ff,
                'neq':{ 
                    'flag': False
                },
                'dt': dt,
                'temp': temp,
                'temp_damp': temp_damp,   # langevin
                'press': press,
                'press_damp': press_damp,
                'neigh_freq': 2,
                'thermo_freq': 2000,
                'style': 'lj/cut/coul/long',
                'ensemble': 'nvt',
                'dump': {
                    'flag': False, 
                    'atom': "all",
                    'freq': 40,
                    'file': "spce_bulk_eq.lammpstrj",
                    'style': "id type xu yu zu #vx vy vz"
                },
                'confinement': confinement.copy(),
                'fixid': ['eqmd','eql'],
                'cutoff':{
                    'lj':   cutoff_lj,
                    'coul': cutoff_coul
                },
                'intramol':{
                    'bond': 'harmonic',
                    'angle': 'harmonic'
                },
                'LJ_mixing' : lj_mixing,
                'LJ_tail'   : lj_tail,
                'LJ_shift'  : lj_shift,
                'shake': {
                    'flag'  :   False,
                    #'fix'   :   "fix fOHs        all    shake 0.0001 20 ${rattle_freq} b 1 a 1",
                    'fix'   :   "fix fOHs        all    shake 0.00001 40 ${rattle_freq} b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
            #########################
            # NEMD
            #########################
            'neqMD':{
                'mol': False,
                'evolve':{
                    'flying'        : "subtract" # example: "substract" for salt-exchange; "union" for water/lj exchange
                },
                'ff': ff,
                ####################
                # H4D
                ####################
                'neq':{
                    "flag": True,               # do NOT change this
                    "max_altitude_neq": 1.0,    # w_max in the paper
                    "max_atom_type_neq": 5,     # for cation; max type of non-flying particles (in this example, 5 types for non-flying)
                    "max_atom_type2_neq": 11,   # for anion; optional
                    "altitude_protocol":{       # only use constant-velocity protocol
                        "nsteps": None,         # will be updated; integer
                        "convex": False,        # do NOT change this
                        "power": 1,             # do NOT change this (constant velocity)
                        # screening parameters
                        "screening_sr": 1.,     # short-range Coulomb screening along altitude (kappa^sr_scr)
                        "screening_lr": 1.      # long-range Coulomb screening along altitude (kappa^lr_scr)
                    },
                    # followings will change along simulations
                    "insertion": None,          # bool; true or false
                    "current_step_neq": None    # integer
                },
                'dt': dt_neq,
                'thermo_freq': 2000, 
                'neigh_freq': 2,
                'style': 'lj/cut/coul/pppm/neq',
                'long_accuracy': 0.0001,
                'ensemble': 'nve',
                'dump': { 
                    'flag': False,
                    'atom': "all",
                    'freq': 600,
                    'file': "spce_bulk_neq.lammpstrj",
                    'style': "id type xu yu zu #vx vy vz"
                },
                'confinement': confinement.copy(),  # NOT used in this version
                'fixid': ['eqmd'],                  # do NOT change this for this example
                'cutoff':{
                    'lj':   cutoff_lj,
                    'coul': cutoff_coul_neq
                },
                'intramol':{
                    'bond': 'harmonic',
                    'angle': 'harmonic'
                },
                'LJ_mixing' :   lj_mixing,
                'LJ_tail'   :   lj_tail,
                'LJ_shift'  :   lj_shift,
                'energies': {           # will be updated during the simulations
                    'ke_new': 0.,
                    'pe_new': 0.,
                    'ke_old': 0.,
                    'pe_old': 0.,
                    'bond_w_old': 0.,
                    'angle_w_old': 0.,
                    'bond_w_new': 0.,
                    'angle_w_new': 0.
                },
                'mom_rev': None,        # will be updated every single gcmc move
                'shake': {
                    'flag'  :   False,
                    'fix'   :   "fix fOHs        all    shake 0.00001 40 0 b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
            'hybrid': {                 # NOT used in this version
                'neq': True,        # non-eq MD trial
                'instant': False,   # instant trial
                'widom': False      # widom insertion
            }
        }

        #############################################################
        # DO NOT NEED TO CHANGE ANYTHING BELOW
        #############################################################

        # update cutoff
        for item in ['lj','coul']:
            self.dic_mc['cutoff'][item]=self.dic_inputs['neqMD']['cutoff'][item]

        # Peratom quantities
        self.dic_peratom={
            "new":{
                "type": None,
                "q": None,
                "x": None,
                "v": None,
                "f": None,
                "image": None
            },
            "old":{
                "type": None,
                "q": None,
                "x": None,
                "v": None,
                "f": None,
                "image": None
            },
            "current":{
                "type": None,
                "q": None,
                "x": None,
                "v": None,
                "f": None,
                "image": None
            }
        }

    def update_noneq_altitude_protocol(self,nstep=100, dic_inputs=None, output="mc_test"):
        self.dic_inputs['neqMD']['neq']["altitude_protocol"]["nsteps"]=nstep

    def update_noneq_velocity(self,nstep=100, dic_inputs=None, output="mc_test"):
        dt=dic_inputs['neqMD']['dt']
        d=dic_inputs['neqMD']['neq']['max_altitude_neq']
        try:
            v=-d/float(nstep)/dt
            print ("= hybrid neMD/MC with v_w=%f"%v,file=output)
        except ZeroDivisionError:
            v=0.
            print ("= Instant addition/deletion",file=output)
            
        dic_inputs['neqMD']['neq']['vspeed_neq']=v

    def update_current_step(self,mylmp=None,dic_inputs=None):
        step_now = mylmp.extract_global("ntimestep",0)
        dic_inputs['neqMD']['neq']['current_step_neq']=step_now
        return step_now

    def check_current_noneq_dic(self,output="dic_test",dic_inputs=None):
        print ("""Test neqMD dictionary
\tcurrent timestep for neqMD: %d
\tchosen move (insertion? %s)"""%(
        dic_inputs['neqMD']['neq']['current_step_neq'],
        dic_inputs['neqMD']['neq']['insertion']
        ), file=output)

    # clear dictionary after trial move
    def clear_dictionary(self,output="dic_test",dic_list=None):
        dic_mc=dic_list["mc"]
        dic_inputs=dic_list["inputs"]
        print ("Before clear energies in noneq dic:",file=output)
        print (dic_inputs['neqMD']['energies'],file=output)
        dic_mc["select_insert"]=None
        dic_mc["select_water"]=None
        dic_mc["selected"]['list']=None
        dic_mc["selected"]['water']=None
        dic_mc["selected"]['salt']['cation']=None
        dic_mc["selected"]['salt']['anion']=None
        dic_mc["accept"]=None
        # clear energies
        for item in dic_inputs['neqMD']['energies'].keys():
            dic_inputs['neqMD']['energies'][item]=0.
        print ("After clear energies in noneq dic:",file=output)
        print (self.dic_inputs['neqMD']['energies'],file=output)

    # replace element
    def replace_dic_peratom(self,r,l,status="current",dic_peratom=None):
        dic_peratom[status]['x'][r]=dic_peratom[status]['x'][l]
        dic_peratom[status]['v'][r]=dic_peratom[status]['v'][l]
        dic_peratom[status]['f'][r]=dic_peratom[status]['f'][l]
        dic_peratom[status]['image'][r]=dic_peratom[status]['image'][l]

    # swap element
    def swap_dic_peratom(self,r,l,status="current",dic_peratom=None):
        # intermediate <- r
        xr=dic_peratom[status]['x'][r]
        vr=dic_peratom[status]['v'][r]
        fr=dic_peratom[status]['f'][r]
        ir=dic_peratom[status]['image'][r]
        # r <- l
        dic_peratom[status]['x'][r]=dic_peratom[status]['x'][l]
        dic_peratom[status]['v'][r]=dic_peratom[status]['v'][l]
        dic_peratom[status]['f'][r]=dic_peratom[status]['f'][l]
        dic_peratom[status]['image'][r]=dic_peratom[status]['image'][l]
        # l <- intermediate
        dic_peratom[status]['x'][l]=xr
        dic_peratom[status]['v'][l]=vr
        dic_peratom[status]['f'][l]=fr
        dic_peratom[status]['image'][l]=ir

    # initialize peratom dictionary
    def init_peratom(self,mylmp=None,dic_peratom=None):
        x,v,f,q,atype,image=mymisc.allocate_xvf(mylmp=mylmp)
        for item in ["new","old","current"]:
            dic_peratom[item]['x']=x
            dic_peratom[item]['v']=v
            dic_peratom[item]['f']=f
            dic_peratom[item]['q']=q
            dic_peratom[item]['type']=atype
            dic_peratom[item]['image']=image

    # update peratom dictionary
    def update_peratom(self,mylmp=None,status="old",dic_peratom=None):
        x,v,f,q,atype,image=mymisc.gather_atoms_xvf(mylmp=mylmp)
        dic_peratom[status]['x']=x
        dic_peratom[status]['v']=v
        dic_peratom[status]['f']=f
        dic_peratom[status]['q']=q
        dic_peratom[status]['type']=atype
        dic_peratom[status]['image']=image

    def update_force(self,mylmp=None,status="current",dic_peratom=None):
        mylmp.scatter_atoms("f",1,3,dic_peratom[status]['f'])            # send it to LAMMPS
        mylmp.command("run 0 # after update force with shake constraint")

    def update_velocity(self,mylmp=None,status="current",dic_peratom=None):
        mylmp.scatter_atoms("x",1,3,dic_peratom[status]['x'])            # send it to LAMMPS
        mylmp.scatter_atoms("v",1,3,dic_peratom[status]['v'])            # send it to LAMMPS
        mylmp.scatter_atoms("image",0,3,dic_peratom[status]['image'])    # send it to LAMMPS
        mylmp.command("run 0 # after update velocity with new forcefield")


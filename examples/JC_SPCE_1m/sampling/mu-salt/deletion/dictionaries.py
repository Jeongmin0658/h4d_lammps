import numpy as np
import math
import os,sys

# source folder
source_where="/SSD/jmkim/Simulation/H4D_method/Github/Github_h4d_Ongoing/ver0.0"
sys.path.append(source_where)

from define_lmp_misc import LMP_misc
mymisc=LMP_misc()

class LMP_dics():
    def __init__(self):
        kJ_to_kcal=1.0/4.184
        ###
        # chemical potential
        # will be updated later
        # ref: Luc Belloni J. Chem. Phys. 149, 094111 (2018)
        # in unit of kT # De Broglie wavelength = 1 A
        chem_pot_water=-15.34   #38.10*kJ_to_kcal
        chem_pot_salt=-314.5    #*thermal # 298.15 K
        # ref: Ind. Eng. Chem. Res. 2017, 56, 4119−4135
        # ref: https://aip.scitation.org/doi/pdf/10.1063/1.453665
        chem_pot_lj=-4.45 
        #chem_pot_salt=-14.85 #-15.94
        # in unit of kT # De Broglie wavelength = 1 A
        # ref:  J. Chem. Phys. 72, 2907 (1980) 
        ###

        # region to insert
        region_to_insert="NULL" # example; "confined" for a certain region, or "NULL" = full box

        # force-field
        ff='JC_SPCE'  # example:'SD_SPCE' # or 'JC_SPCE' or 'simple_lj'
        # unit in lammps
        lj_unit=False # True if lammps uses lj unit

        # temp
        temp=298.15 #K
        press=1.0 #0.2*temp

        ## reduced chemical potential by kT
        ## excess chemical potential to full chemical potential
        #chem_pot_lj=chem_pot_lj/temp+np.log(0.5918) #4.45
        ##chem_pot_lj/=temp

        # timestep
        dt = 2.0 #0.005
        dt_neq = 4.0 #0.002

        # damping
        temp_damp=100.0*dt
        press_damp=1000.0*dt

        # LJ parameter setting
        lj_mixing='arithmetic'
        lj_tail='no'
        lj_shift='yes'

        # cutoff distances
        cutoff_lj = 9.0 #2.5
        cutoff_coul = 9.0 #3.5
        cutoff_lj_neq = 14.0
        cutoff_coul_neq = 14.0

        # implicit confinement
        confinement = {
            'flag': False, #True,
            'wall': "wall/lj93",
            'eps': 8.0,
            'sig': 1.0,
            'cut': 40.0,
            'spacing': 1.0/math.sqrt(2.0)
        }

        # MC inputs
        self.dic_mc={
            "ref"           : "Luc Belloni, J. Chem. Phys. 149, 094111 (2018)",
            "gcmc"          : False, #True,         # if True, no move is accepted, but measure chemical potential
            "loop"          : False,        # not in use
            "lj_flag"       : lj_unit,      # lj units
            "atom_type"     : {             # to keep track of atom list
                "lj"        : [],
                "water"     : [1],          # only oxygen
                "cation"    : [3],
                "anion"     : [4]
            },
            "bias"          :{              # for an ion-pair exchange
                "bias_f"    : "bimodal",
                "alpha"     : 1.0,          # inverse of width
                "center"    : 2.0,          # two center of bimodal
                "ex_vol"    : 0.0,          # excluded volume (length)
                "factor_ins"  : None,       # factor in acceptance prob. with the bias
                "factor_del"  : None,       # factor in acceptance prob. with the bias
                "selection": {              # will be updated if there is a selection process (salt-pair)
                    "ins": 1.,
                    "del": 1.
                }
            },
            # water + lj + salt = 1
            # go down "evolve" in neqMD section, and change it accordingly
            # example: "substract" for salt-exchange; "union" for water/lj exchange
            "prob_water"    : 0.,           # electrolytes
            "prob_lj"       : 0.,           # simple lj solvent
            "prob_salt"     : 1.,           # electrolytes
            "prob_insert"   : 0.0,          # for all ff; either 0.5 or 1 or 0
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
            "confined"      : confinement.copy()["flag"],
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

        # Flying molecule inputs
        self.dic_flying={
            "lj":{ 
                "molecule"  : False, #True,
                "read_input": False, # if True, need to specify how to generate velocity
                "region_to_insert": region_to_insert,
                "natoms"    : 1,
                "group"     : ["ljs"],
                "type_neq_s": "6",    # string
                "type_eq_s" : "1",
                "type_neq_l": [6],   # list
                "type_eq_l" : [1],
                "mass"      : [1.0],  
                "charge"    : [0.],
                "full_list" : { # full list (match with dic_mc[selected][list])
                    "charge": [0.],
                    "type_neq_l": [6],
                    "type_eq_l" : [1]
                },
                "file"      : None, #source_where+"/input_lmp/water_spce_widom.txt",
                "mol_name"  : "widom_lj",
                "activity"  : None #chem_pot_water
            },
            "salt_lj":{ 
                "molecule"  : False, #True,
                "read_input": False, # if True, need to specify how to generate velocity
                "region_to_insert": region_to_insert,
                "natoms"    : 2,
                "group"     : ["ljs"],
                "type_neq_s": "7 8",    # string
                "type_eq_s" : "2 3",
                "type_neq_l": [7,8],   # list
                "type_eq_l" : [2,3],
                "mass"      : [1.0, 1.0],
                "charge"    : [1.0,-1.0],
                "full_list" : { # full list (match with dic_mc[selected][list])
                    "charge": [1.0,-1.0],
                    "type_neq_l": [7,8],
                    "type_eq_l" : [2,3]
                },
                "file"      : None, #source_where+"/input_lmp/water_spce_widom.txt",
                "mol_name"  : "widom_lj",
                "activity"  : None #chem_pot_water
            },
            "water":{
                "molecule"  : True,
                "read_input": True, #False, # if True, need to specify how to generate velocity
                "region_to_insert": region_to_insert,
                "natoms"    : 3,
                "group"     : ["oxygens","hydrogens"],
                "type_neq_s": "8 9",    # string
                "type_eq_s" : "1 2",
                "type_neq_l": [8, 9],   # list
                "type_eq_l" : [1, 2],
                "mass"      : [15.999400, 1.007940, 1.007940],
                "charge"    : [-0.847600,0.423800],
                "full_list" : { # full list (match with dic_mc[selected][list])
                    "charge": [-0.847600,0.423800,0.423800],
                    "type_neq_l": [8,9,9],
                    "type_eq_l" : [1,2,2]
                },
                "file"      : source_where+"/input_lmp/water_spce_widom.txt",
                "mol_name"  : "widom_water",
                "activity"  : None #chem_pot_water
            },
            # salt in use
            "salt":{    # NaCl
                "molecule"  : False,
                "read_input": False,
                "region_to_insert": region_to_insert, # NULL = full box
                "natoms"    : 2,
                "group"     : ["cations","anions"],
                "type_neq_s": "10 11",
                "type_eq_s" : "3 4",
                "type_neq_l": [10, 11],
                "type_eq_l" : [3, 4],
                "mass"      : [22.989769, 35.453000],  
                "charge"    : [1.,-1.],
                "full_list" : { # full list
                    "charge": [1.,-1.],
                    "type_neq_l": [10,11],
                    "type_eq_l" : [3,4]
                },
                "file"      : None,
                "mol_name"  : "widom_salt",
                "activity"  : None #-10
            }
        }

        # update chemical potential
        self.dic_flying['water']['activity']=chem_pot_water
        self.dic_flying['salt']['activity']=chem_pot_salt
        self.dic_flying['lj']['activity']=chem_pot_lj

        # MD input files
        self.dic_inputs={
            'random_seed':  967,
            'input_file':{          # input lammps file to read
                'to_read':  'lmp_jc_re',
                'where':    './'
            },
            'output_file':  'output_test',
            'confinement_factor': 1.0, #3.5/7.5,# for restricted volume
            'cutoff':{
                'lj':       cutoff_lj,
                'coul':     cutoff_coul
            },
            'intramol':{
                'bond':     'harmonic',
                'angle':    'harmonic'
            },
            'LJ_mixing':    lj_mixing,
            'LJ_tail':      lj_tail,
            'LJ_shift':     lj_shift,
            'thermal': None,            # will be updated
            'temp':  temp,
            'temp_damp': temp_damp,       # langevin damping parameter
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
                'ensemble': 'npt',
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
                    'flag'  :   True, #False,
                    'fix'   :   "fix fOHs        all    shake 0.0001 20 ${rattle_freq} b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
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
                'neigh_freq': 4,
                'thermo_freq': 2000,
                'style': 'lj/cut/coul/long',
                'ensemble': 'npt',
                'dump': {
                    'flag': False, 
                    'atom': "all",
                    'freq': 50,
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
                    'flag'  :   True, #False,
                    #'fix'   :   "fix fOHs        all    shake 0.0001 20 ${rattle_freq} b 1 a 1",
                    'fix'   :   "fix fOHs        all    shake 0.00001 40 ${rattle_freq} b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
            'neqMD':{
                'mol': False,
                'evolve':{
                    'flying'        : "subtract" # example: "substract" for salt-exchange; "union" for water/lj exchange
                },
                'ff': ff,
                'neq':{
                    "flag": True,
                    "max_altitude_neq": 3.0,
                    "max_atom_type_neq": 7,
                    "max_atom_type2_neq": 11, # for anion
                    "altitude_protocol":{
                        "nsteps": None,     # will be updated; integer
                        "convex": False,    # True,
                        "power": 1,         # will not change (constant velocity)
                        "screening_sr": 1.,     # coulomb screening along altitude
                        "screening_lr": 1.     # coulomb screening along altitude
                    },
                    # followings will change along simulations
                    "insertion": None,          # bool; true or false
                    "current_step_neq": None    # integer
                },
                'dt': dt_neq,
                'thermo_freq': 4000, #30,
                'neigh_freq': 2,
                'style': 'lj/cut/coul/pppmneq/neq',
                'long_accuracy': 0.0001,
                'ensemble': 'rigid/nve',
                'dump': { 
                    'flag': False,
                    'atom': "all",
                    'freq': 10,
                    'file': "spce_bulk_neq.lammpstrj",
                    'style': "id type xu yu zu #vx vy vz"
                },
                'confinement': confinement.copy(),
                'fixid': ['eqmd','eqions'],
                'cutoff':{
                    'lj':   cutoff_lj_neq,
                    'coul': cutoff_coul_neq
                },
                'intramol':{
                    'bond': 'harmonic',
                    'angle': 'harmonic'
                },
                'LJ_mixing' :   lj_mixing,
                'LJ_tail'   :   lj_tail,
                'LJ_shift'  :   lj_shift,
                'energies': {   # will be updated
                    'ke_new': 0.,
                    'pe_new': 0.,
                    'ke_old': 0.,
                    'pe_old': 0.,
                    # to check water intramolecular energies
                    'bond_w_old': 0.,
                    'angle_w_old': 0.,
                    'bond_w_new': 0.,
                    'angle_w_new': 0.
                },
                'mom_rev': None,    # will be updated every single gcmc move
                'shake': {
                    'flag'  :   False,
                    'fix'   :   "fix fOHs        all    shake 0.00001 40 0 b 1 a 1",
                    'unfix' :   "unfix fOHs"
                }
            },
            'hybrid': { # not in use right now
                'neq': True,        # non-eq MD trial
                'instant': False,   # instant trial
                'widom': False      # widom insertion
            }
        }

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


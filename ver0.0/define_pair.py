from ctypes import *
import sys,random,math
import numpy as np
import os

# dictionary for inputs and data
from dictionaries import LMP_dics
mydics=LMP_dics()
## misc. subroutines
#from define_lmp_misc import LMP_misc
#mymisc=LMP_misc()

# Class for hybrid MD/MC
class LMP_pair():
    def __init__(self):
        self.intro="Define pair potentials"
        # Different altitude for a salt pair
        self.w_add=0. #0.5*(2.583+4.401)*0.5

    # pair potential for equilibrium MD
    def define_force_field(self,mylmp=None,dic=None,dic_peratom=None,loop=False):
        # a non-eq switching protocol
        if dic['neq']['flag']:
            # non-eq dynamics
            step_now=dic['neq']['current_step_neq']
            nsteps=dic['neq']['altitude_protocol']['nsteps']
            max_altitude=dic['neq']['max_altitude_neq']
            max_atom_type=dic['neq']['max_atom_type_neq']
            max_altitude2=dic['neq']['max_altitude_neq']+self.w_add     # for a big anion
            max_atom_type2=dic['neq']['max_atom_type2_neq']             # for a big anion
            lambda_ws=dic['neq']['altitude_protocol']['screening_sr']   # for screening
            lambda_wl=dic['neq']['altitude_protocol']['screening_lr']   # for screening
            if dic['neq']['insertion']:
                isNeqinsertion=1    
            else:
                isNeqinsertion=0    
            # altitude schedule: constant-velocity or convex or concave with a variable exponent
            if dic['neq']['altitude_protocol']['convex']:
                isNeqconvex=1   
            else:
                isNeqconvex=0   
            exponent=dic['neq']['altitude_protocol']['power']
            neq_altitude_protocol="%d %d %f %d %d %d %f %f %d %f"%(step_now,nsteps,max_altitude,max_atom_type,isNeqinsertion,isNeqconvex,exponent,max_altitude2,max_atom_type2,lambda_ws)
            neq_altitude_protocol_long="%d %d %f %d %d %d %f %f %d %f"%(step_now,nsteps,max_altitude,max_atom_type,isNeqinsertion,isNeqconvex,exponent,max_altitude2,max_atom_type2,lambda_wl)
        # define groups only for eq-MD
        if not dic['neq']['flag']: 
            mylmp.commands_string(self.update_groups(ff=dic['ff']))
        ##########################
        # define pair style
        ##########################
        mylmp.command("######################################################################")
        cutlj=dic['cutoff']['lj']
        cutcoul=dic['cutoff']['coul']
        style=dic['style']
        # eqPES
        if style=="lj/cut":
            mylmp.command("pair_style   lj/cut  %f"             % cutlj)
        elif style=="lj/cut/coul/cut":
            mylmp.command("pair_style   lj/cut/coul/cut  %f %f" % (cutlj,cutcoul))
        elif style=="lj/cut/coul/long":
            mylmp.command("pair_style   lj/cut/coul/long %f %f" % (cutlj,cutcoul))
        # time-dependent nePES
        elif style=="lj/cut/neq":
            mylmp.command("pair_style   lj/cut/neq  %f %s" % (cutlj,neq_altitude_protocol))
        elif style=="lj/cut/coul/cut/neq":
            mylmp.command("pair_style   lj/cut/coul/cut/neq  %f %f %s" % (cutlj,cutcoul,neq_altitude_protocol))
        elif style=="lj/cut/coul/pppm/neq" or style=="lj/cut/coul/pppmneq/neq":
            mylmp.command("pair_style   lj/cut/coul/long/neq %f %f %s" % (cutlj,cutcoul,neq_altitude_protocol))
        elif style=="lj/cut/coul/ewald/neq" or style=="lj/cut/coul/ewaldneq/neq":
            mylmp.command("pair_style   lj/cut/coul/long/neq %f %f %s" % (cutlj,cutcoul,neq_altitude_protocol))
        elif style=="lj/cut/coul/dsf/neq":
            mylmp.command("pair_style   lj/cut/coul/dsf/neq 0.2 %f %f %s" % (cutlj,cutcoul,neq_altitude_protocol))
        else:
            return False, "ERORR:::style input is not in the list %s"%style
        # long-range coulomb
        if "long" in style:
            if "neq" in style:
                exit("Long-range Coulomb should be specified either with pppm(neq) or ewald(neq) in NEMD (not **long**)")
            mylmp.command("kspace_style pppm 0.0001")
        elif "pppmneq" in style:
            accuracy=style=dic['long_accuracy']
            mylmp.command("kspace_style pppm/neq %s %s"%(accuracy,neq_altitude_protocol_long))
        elif "ewaldneq" in style:
            mylmp.command("kspace_style ewald/neq 0.0001 %s"%neq_altitude_protocol_long)
        elif "pppm" in style:
            mylmp.command("kspace_style pppm 0.0001")
        elif "ewald" in style:
            mylmp.command("kspace_style ewald 0.0001")
        else:
            mylmp.command("kspace_style none")
        # intramolecular
        if dic['mol']:
            mylmp.commands_string("""
bond_style    %s 
angle_style   %s 
            """%(dic['intramol']['bond'],dic['intramol']['angle']))
        # pair lj parameters
        mylmp.commands_string(self.pair_coeff_setting(ljcut=cutlj, neq=dic['neq']['flag'], ff=dic['ff']))
        # pair modify tail yes or shift yes
        mylmp.command("pair_modify  tail %s shift %s mix %s # LB mixing rule"%(dic['LJ_tail'],dic['LJ_shift'],dic['LJ_mixing']))
        mylmp.command("neigh_modify delay %d every 1 check yes page 500000 one 50000"%dic['neigh_freq'])
        # save internally
        #mylmp.command("run 0 post no # update forcefield (neq? %s)"%dic['neq']['flag'])
        mylmp.command("######################################################################")

    # Groups
    def update_groups(self,ff="JC_SPCE"):
        if ff=='SD_SPCE' or ff=="JC_SPCE":
            return """
group   solvent type 1 2 
group   water   type 1 2 
group   ions    type 3 4
group   wall    type 5 6 7
group   flying  type 8 9 10 11
group   eq_run  union solvent ions flying
            """
        if ff=='simple_lj':
            return """
group   solvent type 1
group   water   empty
group   ions    type 2 3
group   wall    type 4 5
group   flying  type 6 7 8
group   eq_run  union solvent ions flying
            """

    def del_groups(self):
        return """
group   solvent delete
group   water   delete
group   ions    delete
group   wall    delete
group   flying  delete
group   eq_run  delete
        """

    # LJ interaction potential
    def pair_coeff_setting(self,ljcut=9.,neq=False,ff='SD_SPCE'):
        if ff=='JC_SPCE':
            return """
# LJ {ff} (neq? {neq})
# real particles
pair_coeff  1 1     0.155354    3.169   {ljcut}   # O O
pair_coeff  3 3     0.352644    2.1595  {ljcut}   # Na Na
pair_coeff  4 4     0.0127851   4.8305  {ljcut}   # Cl Cl
# flying particles
pair_coeff  8 8     0.155354    3.169   {ljcut}   # O O
pair_coeff  10 10   0.352644    2.1595  {ljcut}   # Na Na
pair_coeff  11 11   0.0127851   4.8305  {ljcut}   # Cl Cl
# wall particles: not exist now
pair_coeff  * 5 0 0 {ljcut}
pair_coeff  * 6 0 0 {ljcut}
pair_coeff  * 7 0 0 {ljcut}
# hydrogens
pair_coeff  * 2 0 0 {ljcut} # hydrogen
pair_coeff  * 9 0 0 {ljcut} # flying hydrogen

# intramolecular force field
bond_coeff  1   100000   1.0
angle_coeff 1   100000   109.470000

        """.format(ljcut=ljcut,neq=neq,ff=ff)
        if ff=='SD_SPCE':
            return """
# LJ {ff} (neq? {neq})
# real particles
pair_coeff  1 1     0.155354    3.169   {ljcut}   # O O
pair_coeff  3 3     0.100039    2.583   {ljcut}   # Na Na
pair_coeff  4 4     0.100039    4.401   {ljcut}   # Cl Cl
# flying particles
pair_coeff  8 8     0.155354    3.169   {ljcut}   # O O
pair_coeff  10 10   0.100039    2.583   {ljcut}   # Na Na
pair_coeff  11 11   0.100039    4.401   {ljcut}   # Cl Cl
# wall particles: not exist now
pair_coeff  * 5 0 0 {ljcut}
pair_coeff  * 6 0 0 {ljcut}
pair_coeff  * 7 0 0 {ljcut}
# hydrogens
pair_coeff  * 2 0 0 {ljcut} # hydrogen
pair_coeff  * 9 0 0 {ljcut} # flying hydrogen

# intramolecular force field
bond_coeff  1   100000   1.0
angle_coeff 1   100000   109.470000

        """.format(ljcut=ljcut,neq=neq,ff=ff)
        if ff=='simple_lj':
            return """
# LJ (neq? {neq})
# real lj particles
pair_coeff  * *         1.  1.  {ljcut}   # solvent + ions
pair_coeff  4*5 4*5     0   0   {ljcut}   # wall atoms

        """.format(ljcut=ljcut,neq=neq)

    # for delete trial move: flying particles are non-interacting
    def pair_coeff_delete_trial(self,ljcut=9.,ff='SD_SPCE'):
        if ff=='JC_SPCE':
            return """
# Non-interacting flying molecule (trial deletion)
# real particles
pair_coeff  1 1     0.155354    3.169   {ljcut}   # O O
pair_coeff  3 3     0.352644    2.1595   {ljcut}   # Na Na
pair_coeff  4 4     0.0127851   4.8305   {ljcut}   # Cl Cl
# flying particles
pair_coeff  * 8  0 0 {ljcut}    # non-interacting
pair_coeff  * 10 0 0 {ljcut}    # non-interacting
pair_coeff  * 11 0 0 {ljcut}    # non-interacting
# wall particles: not exist now
pair_coeff  * 5 0 0 {ljcut}
pair_coeff  * 6 0 0 {ljcut}
pair_coeff  * 7 0 0 {ljcut}
# hydrogens
pair_coeff  * 2 0 0 {ljcut}     # hydrogen
pair_coeff  * 9 0 0 {ljcut}     # flying hydrogen

# intramolecular force field
bond_coeff  1   100000   1.0
angle_coeff 1   100000   109.470000

# flying water or ions
set type 8  charge 0.           # non-interacting
set type 9  charge 0.           # non-interacting
set type 10 charge 0.           # non-interacting
set type 11 charge 0.           # non-interacting
        """.format(ljcut=ljcut)
        if ff=='SD_SPCE':
            return """
# Non-interacting flying molecule (trial deletion)
# real particles
pair_coeff  1 1 0.155354    3.169   {ljcut}   # O O
pair_coeff  3 3 0.100039    2.583   {ljcut}   # Na Na
pair_coeff  4 4 0.100039    4.401   {ljcut}   # Cl Cl
# flying particles
pair_coeff  * 8  0 0 {ljcut}    # non-interacting
pair_coeff  * 10 0 0 {ljcut}    # non-interacting
pair_coeff  * 11 0 0 {ljcut}    # non-interacting
# wall particles: not exist now
pair_coeff  * 5 0 0 {ljcut}
pair_coeff  * 6 0 0 {ljcut}
pair_coeff  * 7 0 0 {ljcut}
# hydrogens
pair_coeff  * 2 0 0 {ljcut}     # hydrogen
pair_coeff  * 9 0 0 {ljcut}     # flying hydrogen

# intramolecular force field
bond_coeff  1   100000   1.0
angle_coeff 1   100000   109.470000

# flying water or ions
set type 8  charge 0.           # non-interacting
set type 9  charge 0.           # non-interacting
set type 10 charge 0.           # non-interacting
set type 11 charge 0.           # non-interacting
        """.format(ljcut=ljcut)
        if ff=='simple_lj':
            return """
# Non-interacting flying lj molecule (trial deletion)
# real lj particles
pair_coeff  * *         1.  1.  {ljcut}   # solvent + ions
pair_coeff  4*5 4*5     0   0   {ljcut}   # wall atoms
# flying particles
pair_coeff  * 6         0   0   {ljcut}   # flying solvent  non-interacting
pair_coeff  * 7         0   0   {ljcut}   # flying cation   non-interacting
pair_coeff  * 8         0   0   {ljcut}   # flying anion    non-interacting

# flying solvetn or ions
set type 6  charge 0.                       # non-interacting
set type 7  charge 0.                       # non-interacting
set type 8  charge 0.                       # non-interacting
        """.format(ljcut=ljcut)

    def confinement_wall(self,group="all",on=True,neq=False,wall="wall/lj93",eps=8.0,sig=1.0,cut=15.0,spacing=1.0):
        if on:
            if wall=="wall/lj93":
                m="""
# wall potential
fix walllo {group} {wall} zlo -${box_size_z} {eps} {sig} {cut} units box
fix wallhi {group} {wall} zhi  ${box_size_z} {eps} {sig} {cut} units box
                """.format(group=group,wall=wall,eps=eps,sig=sig,cut=cut,box_size_z="{box_size_z}")
            if wall=="wall/lj1043":
                m="""
# wall potential
fix walllo {group} {wall} zlo -${box_size_z} {eps} {sig} {cut} {spacing} units box
fix wallhi {group} {wall} zhi  ${box_size_z} {eps} {sig} {cut} {spacing} units box
                """.format(group=group,wall=wall,eps=eps,sig=sig,cut=cut,spacing=spacing,box_size_z="{box_size_z}")
            if not neq:
                return m+"""
# turn this on only for eq PES {neq}
fix_modify walllo energy yes 
fix_modify wallhi energy yes 
                """.format(neq=neq)
            else:
                return m
        else:
            return """
unfix walllo
unfix wallhi
            """
    # LJ tail correction
    def LJ_mu_tail(self,rho=1.):
        eps,sig,cut=self.pair[0],self.pair[1],self.pair[2]
        rc3inv=np.power(sig/cut,3.)
        rc9inv=np.power(rc3inv,3.)
        utail=4.*8./3.*math.pi*rho*eps*np.power(sig,3.)*(2.*rc9inv/3. - rc3inv)
        return utail

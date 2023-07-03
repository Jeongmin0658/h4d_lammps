from ctypes import *
import sys,random,math
import numpy as np
import os

## Dictionary as input files
#from dictionaries import LMP_dics
#mydics=LMP_dics()

class LMP_flying():
    def __init__(self):
        self.intro="Define flying particles"

    def set_list_flying_molecules(self,mylmp=None, region_to_insert="confined", mymove=None, dic_list=None, output=None):
        dic_mc=dic_list["mc"]
        dic_flying=dic_list["flying"]

        # update list
        mymove.update_list_atoms(mylmp=mylmp,dic_mc=dic_mc,output=output)

        # define molecules to read
        for flying in ['water','salt','lj']:
            if dic_flying[flying]['read_input'] and dic_mc["prob_"+flying]>0.:
                widom="widom_"+flying
                mylmp.command("molecule     %s    %s"%(widom,dic_flying[flying]['file']))
        #mylmp.command("molecule      widom_salt     %s"%dic_flying['salt']['file'])

        # region to insert
        dic_mc['region_to_insert']=region_to_insert

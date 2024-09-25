import sys
import string
import math
import constants as cnst
from properties import *

class Ref_Rel_E(object):
    def __init__(self, mol_list, ref_data): 
        self.mol_list = mol_list
        self.ref_data = ref_data
        self.ref_reles = []      # list of lists, contains info for all displacements
        self.eq_geos = []        # list of lists, contains info for all equilibrium geometries
        self.gen_ref_rel()          

    def gen_ref_rel(self):
        self.eq_geos = []
        for i in range(len(self.mol_list)):
            if 'eq' in self.mol_list[i]:
	       eq_data = []
               eq_data.append(self.mol_list[i])                           # index 0: filename of this eq geo  
               eq_data.append(self.ref_data.ref_list[i].hof[0].ref_heat)  # index 1: reference heat of this eq geo
               eq_data.append(i)                                          # index 2: the index in mol_list for this eq_geo                
               self.eq_geos.append(eq_data)        

        self.ref_reles = []    
        for i in range(len(self.mol_list)):
            if '01_eq' not in self.mol_list[i]:   #string subject to change
                for j in range(len(self.eq_geos)):
                    if self.eq_geos[j][0] == self.ref_data.ref_list[i].hof[0].eq_geo:
                         rele_data = [] 
                         ref_rele = self.ref_data.ref_list[i].hof[0].ref_heat - self.eq_geos[j][1]
                         rele_data.append(self.mol_list[i])       # index 0: filename for this disp
                         rele_data.append(ref_rele)               # index 1: reference relative energy
                         rele_data.append(i)                      # index 2: the index in mol_list for this disp
                         rele_data.append(self.eq_geos[j][0])     # index 3: the eq geo of this disp
                         rele_data.append(self.ref_data.ref_list[i].hof[0].rela_weight)   # index 4: the weight of this data point
                         self.ref_reles.append(rele_data) 


class PM_Rel_E(object):
    def __init__(self, out_heat, ref_rel):
        self.pm_heat = out_heat
        self.ref_rel = ref_rel
        self.pm_reles = []     # stores semiempirical relative energy
        self.rele_error = []
        self.gen_pm_rel()

    def gen_pm_rel(self):
        self.pm_reles = []
        for i in range(len(self.ref_rel.ref_reles)):
            eq_geo = self.ref_rel.ref_reles[i][3]
            heat = self.pm_heat[self.ref_rel.ref_reles[i][2]]
            if heat != 0:    
                for j in range(len(self.ref_rel.eq_geos)):
                    if eq_geo == self.ref_rel.eq_geos[j][0]:
                         pm_rele = heat - self.pm_heat[self.ref_rel.eq_geos[j][2]]
                         self.pm_reles.append(pm_rele)
                         self.rele_error.append(pm_rele - self.ref_rel.ref_reles[i][1]) 
            else:
                pm_rele = 0
                self.pm_reles.append(pm_rele)
                penalty = 10000
                self.rele_error.append(penalty)


 

import sys
import string
import math
import constants as cnst
from properties import *

###############################################################
# Class to handle reference data and comparisons between      #
# computed data and reference data.                           #
#                                                             #
# Class Ref_Set: (added 4/9/2020)                             #
#   Contains ref_list with list of reference data for all     #
#        molecules in the reference list                      #
#   1. __init__: Store list of molecules & call read_refs     #
#   2. read_refs: For each molecule, set up Ref object to     #
#        store in ref_list                                    #
#                                                             #
# Class Ref:                                                  #
#   Contains all reference data (only excited states for now, #
#        may expand to include other data types later         #
#   1. __init__: Accepts the input line with reference data   #
#        and calls parsing                                    #
#   2. parse: Splits the input line into parts and sets up    #
#        Ex_Ref items for each excited state                  #
#   3. gen_types: Generate a list of AO types needed to read  #
#        in the orbitals for this system                      #
#                                                             #
# Class Ex_Ref:                                               #
#   Contains reference data for one excited state             #
#   1. __init__: Sets up a reference excited state            #
#   2. match_ex_ref: Finds the relevant excited state for     #
#        comparison to the ref data                           #
#   3. print_ex_ref: Prints the reference data, computed      #
#        data, and error                                      #
###############################################################

class Ref_Set(object):
    def __init__(self, mol_list):
        self.mol_list = mol_list
        self.ref_list = []

        self.read_refs()

    def read_refs(self):
        for mol in self.mol_list:
            #mol = mol[:-1]
           # if '.' in mol:
               # mol = mol[:mol.rfind('.')]
        
            # Read .dat file and extract reference data
            dat = open(mol + ".dat",'r')
            dline = dat.readline()
            dline = dat.readline()
            dline = dat.readline()
            refs = Ref(dline)
            self.ref_list.append(refs)


class Ref(object):
    def __init__(self, refline):
        self.refline = refline

        self.exc = []
        self.types = []
        # Add heat of formation as reference data
        self.hof = []

        self.parse()
        self.gen_types()

    def parse(self):
        line = self.refline
        self.exc=[]
        
        # Get heat of formation
        self.hof = []
        if 'H=' in self.refline:
            heat = float(line[line.find('H=')+2:].split()[0])            
            hwt = 1.0  # set a defualt for all
            if 'EG=' in self.refline:
                eqgeo = line[line.find('EG=')+3:].split()[0]
                rwt = float(line[line.find('RW=')+3:].split()[0])
                self.hof.append(Hof_Ref(heat, eqgeo, hwt, rwt))
            elif 'RG=' in self.refline:
                eqgeo = line[line.find('RG=')+3:].split()[0]
                rwt = float(line[line.find('RW=')+3:].split()[0])
                self.hof.append(Hof_Ref(heat, eqgeo, hwt, rwt))

        # Look for up to 9 excited states
        for i in range(1,10):
            base = "E"+str(i)
            if (base + "=") in self.refline:
                #print 'ref ',i
                # Set up default values
                dip    = ''
                dipwt  = 0.6
                occ    = ''
                occwt  = 0.6
                vir    = ''
                virwt  = 0.6
                weight = 1.0
                source = 'none'

                energy = float(line[line.find(base + "=")+3:].split()[0])
                if (base + "d=") in self.refline:
                    dip    = line[line.find(base + "d=")+4:].split()[0]
                    if (base + "dw=") in self.refline:
                       dipwt  = line[line.find(base + "dw=")+4:].split()[0]
                if (base + "o=") in self.refline:
                    occ    = line[line.find(base + "o=")+4:].split()[0]
                    if (base + "ow=") in self.refline:
                       occwt  = float(line[line.find(base + "ow=")+5:].split()[0])
                if (base + "v=") in self.refline:
                    vir    = line[line.find(base + "v=")+4:].split()[0]
                    if (base + "vw=") in self.refline:
                       virwt  = float(line[line.find(base + "vw=")+5:].split()[0])
                if (base + "w=") in self.refline:
                    weight = float(line[line.find(base + "w=")+4:].split()[0])
                if (base + "s=") in self.refline:
                    source = line[line.find(base + "s=")+4:].split()[0]

                #print energy, dip, dipwt, occ, occwt, vir, virwt
                #print line

                self.exc.append(Ex_Ref(energy, dip, dipwt, occ, occwt, vir, virwt, weight, source))

    def gen_types(self):
        self.types = []

        # Add a types line for occ & vir for every reference to simplify later counting
        for ref in self.exc:
            if len(ref.occ) > 0:
                occ = ref.occ.split('_')
                self.types.append(occ)
                if len(occ) > 1:
                    self.types[-1].append('Orbtype')
                elif 'S' in ref.occ or 'P' in ref.occ or 'D' in ref.occ:
                    self.types[-1].append('Orb')
                else:
                    self.types[-1].append('Attype')
            else:
                self.types.append([])

            if len(ref.vir) > 0:
                vir = ref.vir.split('_')
                self.types.append(vir)
                if len(vir) > 1:
                    self.types[-1].append('Orbtype')
                elif 'S' in ref.occ or 'P' in ref.occ or 'D' in ref.occ:
                    self.types[-1].append('Orb')
                else:
                    self.types[-1].append('Attype')
            else:
                self.types.append([])
        #print self.types
            
    def match_refs(self,states):
        for i in range(0,len(self.exc)):
            ref = self.exc[i]
            ref.match_ex_ref(states,i)

    def calc_error(self,pm_heat):
        for i in range(0,len(self.hof)): 
            ref = self.hof[i]
            ref.error_hof(pm_heat)


class Ex_Ref(object):
    def __init__(self, energy, dip, dipwt, occ, occwt, vir, virwt, weight, source):
        self.energy = energy
        self.dip    = dip
        self.dipwt  = dipwt
        self.occ    = occ
        self.occwt  = occwt
        self.vir    = vir
        self.virwt  = virwt
        self.weight = weight
        self.source = source

    def match_ex_ref(self,states,refnum):
        # states contains excited-state data read in from the Mopac output file
        # Loop over states from lowest to highest, starting at 1 (0 is ground state)
        self.match_found = False
        self.match = 0
        i = 1
        dip_cutoff = self.dipwt
        dip_zero_cutoff = self.dipwt
        occ_cutoff = self.occwt
        vir_cutoff = self.virwt

        while self.match_found == False and i < len(states):
            curr = states[i]
            dip_match = False
            occ_match = False
            vir_match = False
            if curr.ref_matched == False:
                if len(self.dip) == 0:
                    dip_match = True
                    #print i, 'No dip ref'
                else:
                    tmag = math.sqrt(curr.tdip[0]**2 + curr.tdip[1]**2 + curr.tdip[2]**2) + 0.0001
                    #print i, curr.tdip, tmag, self.dip.lower()
                    if   'x' in self.dip.lower() and abs(curr.tdip[0]/tmag) > dip_cutoff:
                         dip_match = True
                    elif 'y' in self.dip.lower() and abs(curr.tdip[1]/tmag) > dip_cutoff:
                         dip_match = True
                    elif 'z' in self.dip.lower() and abs(curr.tdip[2]/tmag) > dip_cutoff:
                         dip_match = True
                    elif '0' in self.dip.lower() and tmag < dip_zero_cutoff:
                         dip_match = True
                    #print i, curr.tdip, tmag, self.dip.lower(), dip_match

                if len(self.occ) == 0:
                    occ_match = True
                    #print i, 'No occ ref'
                else:
                    occ_type = refnum * 2
                    if curr.ochar[occ_type] > occ_cutoff:
                         occ_match = True
                    #print i, occ_type, curr.ochar[occ_type], occ_cutoff, occ_match

                if len(self.vir) == 0:
                    vir_match = True
                    #print i, 'No vir ref'
                else:
                    vir_type = occ_type + 1
                    if curr.vchar[vir_type] > vir_cutoff:
                         vir_match = True
                    #print i, vir_type, curr.vchar[vir_type], vir_cutoff, vir_match

                if dip_match and occ_match and vir_match:
                    self.match = i
                    self.match_found = True
                    curr.ref_matched = True
                    #print self.match, refnum
                    #print 'Energy ',self.energy, curr.energy
                #print i, dip_match, occ_match, vir_match, curr.tdip, self.dip.lower(), curr.ochar[occ_type], occ_cutoff, curr.vchar[vir_type], vir_cutoff
            i = i + 1
        self.comp_energy = curr.energy
        self.error_energy = curr.energy - self.energy
        if self.match_found == False:
            #print "Error: No match found for ref " + str(refnum+1)
            self.comp_energy = 0.0
            self.error_energy = 0.0
        #self.print_ex_ref()

    def print_ex_ref(self):
        print 'Energy ', string.rjust('%.3f'%self.energy,8), string.rjust('%.3f'%self.comp_energy,8), string.rjust('%.3f'%self.error_energy,8), string.rjust('%.3f'%self.weight,8)

class Hof_Ref(object):
    def __init__(self, heat, eqgeo, hwt, rwt):
        self.ref_heat = heat
        self.heat_weight = hwt
        self.rela_weight = rwt
        self.eq_geo = eqgeo  
 
    def error_hof(self, pm_heat):
        self.comp_heat = pm_heat
        if self.comp_heat != 0:
            self.match_found = True
            self.error_heat = self.comp_heat - self.ref_heat
        else:
            self.match_found = False
            self.comp_heat = 0.0
            self.error_heat = 0.0



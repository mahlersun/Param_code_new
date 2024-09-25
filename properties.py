import math
import constants as cnst

#############################
#
#  Classes to store molecular properties
#  
#  Classes included:
#     Atom
#     MolOrb
#     Config
#     State 
#
#  RLG, 4/25/17
#
#############################

#### Atom
# coord = [x,y,z]
# elem = string of element abbreviation
# ao = list of AOs with a string of True/False for type information
# mass = atomic mass in a.u. from dictionary
# aotype = list of AO types (S, Px, ...)
####
class Atom(object):
    def __init__(self, coord, elem):
        self.coord = coord
        self.elem = elem
        self.ao = []
        self.mass = cnst.mass_dict[self.elem] 
        self.aotype = []

    # Add atomic orbitals
    def add_ao(self, type):
        self.ao.append(type)

    def add_aotype(self, aotype):
        self.aotype.append(aotype)

#### Molecular orbital
# energy = MO eigenvalue in eV
# coeff = list of AO coefficients
# occ = occupation number (not used??)
# char = list of MO character corresponding to each type in types file 
#        Ranges from 0 to 1, one char element per type
#        Example: Ag SP character
####
class MolOrb(object):
    def __init__(self, energy=0.0, occ=0.0):
        self.energy = energy
        self.coeff = []
        self.occ = occ
        self.char = []

    # Add more AO coefficients
    def add_coeff(self, new_coeff):
        self.coeff.append(float(new_coeff))

#### Electron configuration
# occ = list of occupied orbitals electrons are excited from
# vir = list of virtual orbitals electrons are excited to
# energy = configuration energy in eV
# char = change in character upon excitation corresponding to to each type in types file
#        Ranges from 0 to 1, one char element per type
#        Example: Change in Ag SP character upon excitation
#
# Properties defined elsewhere and rarely used:
#   achar = local character of AOs within the type
#           Example: Fraction of excitation that is Ag SP -> Ag SP
#   bchar = local character of AOs outside of the type
#           Example: Fraction of excitation both from and to non-AgSP orbitals
####
class Config(object):
    def __init__(self, occ, vir, char, energy=0.0):
        self.change_occ(occ)
        self.change_vir(vir)
        self.char = char
        self.energy = energy

    # Set energy of the configuration
    def change_energy(self,energy):
        self.energy = energy

    # Change occupied orbital(s) this configuration excites from
    def change_occ(self, occ):
        if occ is list:
            self.occ = occ
        else:
            self.occ = [occ]

    # Change virtual orbital(s) this configuration excites to
    def change_vir(self, vir):
        if vir is list:
            self.vir = vir
        else:
            self.vir = [vir]

#### Excited state
# energy = excited-state energy in eV (ground state is 0.0)
# osc = oscillator strength (unitless)
# tdip = transition dipole moment (Debye)
# coeff = list of [config_number,config_contribution] 
#         Includes percent contributions (squares of CI coefficients)
# char = change in character upon excitation corresponding to to each type in types file
#        Ranges from 0 to 1, one char element per type
#        Example: Change in Ag SP character upon excitation
#
# Properties defined elsewhere and rarely used:
#   achar = local character of AOs within the type
#           Example: Fraction of excitation that is Ag SP -> Ag SP
#   bchar = local character of AOs outside of the type
#           Example: Fraction of excitation both from and to non-AgSP orbitals
# 
#   energy_shift = difference between the weighted average contributing configuration energy and the actual excited-state energy
#   weighted_std = weighted standard deviation of the contributing configuration energies
#   superatom_char = superatomic character of the orbitals
####
class State(object):
    def __init__(self, energy=0.0, osc=0.0):
        self.energy = energy
        self.osc = osc
        self.coeff = []
        self.char = []
        self.tdip = [0.0,0.0,0.0]

        self.energy_shift = 0.0
        self.weighted_std = 0.0
        self.superatom_char = 0.0

        # For Param - track which excited states have been matched to reference data
        self.ref_matched = False

    def add_coeff(self, config, coeff):
        self.coeff.append([config, coeff])

    def calc_tdip(self, pol):
        for i in range(0,3):
            self.tdip[i] = pol[i]*math.sqrt(self.osc/(cnst.fosc_fact*self.energy))

    def calc_osc(self, tdip):
        self.tdip = tdip
        self.osc = (tdip[0]**2 + tdip[1]**2 + tdip[2]**2) * (cnst.fosc_fact*self.energy)
        #print self.energy, self.tdip

    # Find energy shift related properties
    # Based on Adam's PEW script
    def find_energy_shift(self, configs, mos):
        energies = []
        orb_energies = []
        conf_weights = []
        for conf in self.coeff:
            #print conf
            energies.append(configs[conf[0]-1].energy)
            orb_energies.append((mos[configs[conf[0]-1].occ[0]].energy - mos[configs[conf[0]-1].vir[0]].energy))
            conf_weights.append(conf[1])
            #print conf[0],configs[conf[0]-1].occ[0],configs[conf[0]-1].vir[0],mos[configs[conf[0]-1].occ[0]].energy,mos[configs[conf[0]-1].vir[0]].energy

        # Compute weighted average energy
        # For direct comparison with Adam, base on orbital energies, not configuration energies
        weighted_avg_energy = 0.0
        wt_sum = 0.0
        for i in range(0,len(energies)):
            weighted_avg_energy += orb_energies[i]*conf_weights[i]
            wt_sum += conf_weights[i]
        weighted_avg_energy /= wt_sum
        #self.weighted_average_energy = numpy.average(energies, weights=conf_weights)

        # Compute weighted standard deviation
        # For direct comparison with Adam, base on orbital energies, not configuration energies
        self.weighted_std = 0.0
        for i in range(0,len(orb_energies)):
            self.weighted_std += conf_weights[i] * (orb_energies[i] - weighted_avg_energy)**2
        self.weighted_std /= wt_sum
        self.weighted_std = math.sqrt(self.weighted_std)
        #self.weighted_std = numpy.average((energies - self.weighted_average_energy)**2, weights=conf_weights)

        self.energy_shift = self.energy - weighted_avg_energy

    # Find the collectivity index
    # Based on Adam's PEW script
    def find_collectivity(self):
        participation_sum = 0
        for conf in self.coeff:
            participation_sum += ((conf[1])**2)
        tipr = 1 / participation_sum
        self.collectivity = tipr



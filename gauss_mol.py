import sys
import string
import math
import constants as cnst
from properties import *
from molecule import *

# Holder for all properties of a system
class GaussMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('Charge')
        self.charge = int(self.line.split()[2])

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"


    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read or self.charge_read:
            self.file.seek(0)

        # Check whether the output is an optimization (read optimized geometry) or not (read initial geometry)
        self.read_to('Gaussian 16:')
        self.read_to('---------------')
        self.read_for(1)
        opt = True if 'opt' in self.line.lower() else False

        self.atoms = []
        self.at_types = []

        if opt:
            self.read_to("Optimized Parameters")
        self.read_to("Input orientation")
        self.read_for(5)
        while '-----------' not in self.line and len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[3]),float(line[4]),float(line[5])],self.atomic_symbol.keys()[self.atomic_symbol.values().index(int(line[1]))]))
            if self.atomic_symbol.keys()[self.atomic_symbol.values().index(int(line[1]))] not in self.at_types:
                self.at_types.append(self.atomic_symbol.keys()[self.atomic_symbol.values().index(int(line[1]))])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
            self.atoms_read = True
            #print "Atoms read"
        else:
            print "Error: Atoms not read"

    # Read the number of molecular orbitals
    def read_norbs(self):
        self.numorb = 0
        print 'Not yet implemented for Qchem output'

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        if self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        print 'Not yet implemented for Qchem output'

        # Save that dipole has been read
        if 'Dipole moment=' in self.line:
            self.dipole_read = True
            #print "Dipole read"
        else:
            print "Error: Ground-state dipole not read"

    # Read the molecular orbitals
    def read_orbs(self, types=[], mo_min=1, mo_max=0):
        # Make sure atoms and norbs have been read
        if not self.atoms_read:
            self.read_atoms()
        if not self.norbs_read:
            self.read_norbs()
        if self.configs_read or self.states_read:
            self.file.seek(0)

        print 'Not yet implemented for Qchem output'
        orb_count = 0
        self.mos = []

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print "Orbitals read"
        else:
            print "Error: Orbitals not read"


    # Read electron configurations for CI
    def read_configs(self):
        if self.states_read:
            self.file.seek(0)

        if not self.orbs_read:
            self.read_orbs()

        print 'Not yet implemented for Qchem output'
        self.configs = []

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"


    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
        if not self.configs_read:
            self.read_configs()

        print 'Not yet implemented for Qchem output'
        self.states = [State(0.0,0.0)]

        # Save that states have been read
        if len(self.states) > 0:
            self.states_read = True
            #print "Excited states read"
        else:
            print "Error: Excited states not read"

    def write_vibs(self):
        if not self.atoms_read:
            self.read_atoms()

        self.modes = []
        nmodes = 0
        
        nmo  = open(     'nmodes.inp', 'w')
        while 'STANDARD THERMODYNAMIC' not in self.line and len(self.line) > 0:
            self.read_to('Frequency:')
            nmo.write('            ' + self.line[12:] + '           --------------------   --------------------   --------------------\n')
            self.read_for(7)
            while 'TransDip' not in self.line and len(self.line) > 0:
                nmo.write(self.line)
                self.line = self.file.readline()
            self.read_for(2)
            nmo.write('\n\n')
        nmo.close()

        atm  = open( 'atomicmass.inp', 'w')
        for i in range(0,len(self.atoms)):
            atm.write(string.rjust(str(i),5) + '. ' + string.ljust(self.atoms[i].elem,3) + string.rjust('%.8f'%self.atoms[i].mass,12) + '\n')
        atm.close()





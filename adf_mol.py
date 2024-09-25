import sys
import string
import math
import constants as cnst
from properties   import *
from molecule     import *

# Holder for all properties of a system
class AdfMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        
        self.read_to('title ')
        self.read_to('CHARGE')
        self.charge = int(self.line.split()[1])

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        if self.charge_read == False:
            self.read_to('title input')
        
        # Check if this is a geometry optimization; if so, read optimized geometry
        self.read_to('GEOMETRY')
        opt = False
        while 'END' not in self.line:
            self.line = self.line.lower()
            if 'optim' in self.line:
                opt = True
            self.line = self.file.readline()

        self.atoms = []
        self.at_types = []

        if opt:
            self.read_to('GEOMETRY CONVERGED')
            self.read_for(4)
            while '<' not in self.line:
                line = self.line.split()
                if len(self.atoms) < 9:
                    elem = line[0][2:]
                else:
                    elem = line[0][3:]
                self.atoms.append(Atom([float(line[1]),float(line[2]),float(line[3])],elem))
                if elem not in self.at_types:
                    self.at_types.append(elem)
                self.line = self.file.readline()

            # Reset position in file
            self.file.seek(0)
            self.read_to('title ')

        else:
            self.read_to('ATOMS')
            self.line = self.file.readline()
            
            while 'END' not in self.line and len(self.line) > 0:
                line = self.line.split()
                self.atoms.append(Atom([float(line[1]),float(line[2]),float(line[3])],line[0]))
                if line[0] not in self.at_types:
                    self.at_types.append(line[0])
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
        print "Not implemented for ADF"

        # Save that norbs has been read
        if self.numorb > 0:
            self.norbs_read = True
            #print "Orbital count read"
        else:
            print "Error: Number of orbitals not read"

    # Read the ground state dipole moment into gdip
    def read_dipole(self):
        print "Not implemented for ADF"

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

        self.mos = []
        orb_count = 0
        print "Not implemented for ADF"

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print "Orbitals read"
        else:
            print "Error: Orbitals not read"


    # Read electron configurations for CI
    def read_configs(self,types):
        if self.states_read:
            self.file.seek(0)

        self.configs = []
        print "Not implemented for ADF"

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"


    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
        #if not self.configs_read:
        #    self.read_configs()

        self.read_to('Excitations calculated with the davidson method')
        self.line = self.file.readline()
        nirrep = 0
        while '                                                               ' not in self.line and len(self.line) > 0:
            nirrep += 1
            self.line = self.file.readline()
        print 'nirrep', nirrep

        self.states = [State(0.0,0.0)]

        # Read states of each symmetry
        for i in range(0,nirrep):
            istart = len(self.states)
            self.read_to('Excitation energies E')
            self.read_for(5)
            allowed = False

            while len(self.line) > 2:
                line = self.line.split()
                e = float(line[2])
                o = float(line[3])
                if o > 0.0:
                    allowed = True
                self.states.append(State(e,o))
                self.line = self.file.readline()

            if allowed == True:
                self.read_to('Transition dipole moments')
                self.read_for(5)
                while len(self.line) > 2:
                    line = self.line.split()
                    self.states[istart].tdip = [float(line[3])*4.8032,float(line[4])*4.8032,float(line[5])*4.8032]
                    istart += 1
                    self.line = self.file.readline()


        # Sort states by energy
        self.states.sort(key = lambda x: x.energy)

        #for i in range(0,11):
        #    print self.states[i].energy, self.states[i].tdip

        # Save that states have been read
        if len(self.states) > 1:
            self.states_read = True
            print len(self.states), "excited states read"
        else:
            print "Error: Excited states not read"

    def write_vibs(self):
        if not self.atoms_read:
            self.read_atoms()

        self.modes = []
        nmodes = 0

        nmo  = open(     'nmodes.inp', 'w')
        self.read_to('Vibrations and Normal Modes')
        self.read_for(7)

        while 'List of All Frequencies' not in self.line and len(self.line) > 0:
            nmo.write(self.line)
            self.line = self.file.readline()
        nmo.close()

        atm  = open( 'atomicmass.inp', 'w')
        for i in range(0,len(self.atoms)):
            atm.write(string.rjust(str(i),5) + '. ' + string.ljust(self.atoms[i].elem,3) + string.rjust('%.8f'%self.atoms[i].mass,12) + '\n')
        atm.close()



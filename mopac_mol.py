import sys
import string
import math
import constants as cnst
from properties import *
from molecule import *

# Holder for all properties of a system
class MopacMolecule(Molecule):
    # Read the molecular charge
    def read_charge(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
        
        self.read_to('MOPAC   7.10') # modify this 
        self.read_to('CHARGE','**************************************') # modify this
        self.charge = int(self.line.split()[-1]) if 'CHARGE' in self.line else 0

        # Save that atoms have been read
        self.charge_read = True
        #print "Charge read"

    # Read the atomic coordinates and set up atoms
    def read_atoms(self):
        if self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('CARTESIAN COORDINATES')
        self.read_for(4)
        self.atoms = []
        self.at_types = []
        
        while len(self.line) > 1:
            line = self.line.split()
            self.atoms.append(Atom([float(line[2]),float(line[3]),float(line[4])],line[1]))
            if line[1] not in self.at_types:
                self.at_types.append(line[1])
            self.line = self.file.readline()

        # Save that atoms have been read
        if len(self.atoms) > 0:
            self.atoms_read = True
            #print "Atoms read"
        else:
            print "Error: Atoms not read"

    # Read the number of molecular orbitals
    def read_norbs(self):
        if self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)

        self.read_to('SHELL')
        self.read_for(1)
        line = self.line.split()
        self.homo = 0
        for i in range(1,len(line)):
            self.homo += int(line[i])


        self.read_to('TOTAL')
        line = self.line.split()
        self.numorb = 0
        for i in range(1,len(line)):
            self.numorb += int(line[i])

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

        self.read_to('Dipole moment=')
        self.gdip = float(self.line.split()[-2])
        #print self.gdip

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


        if mo_max == 0:
            mo_max = self.numorb
        self.mos = []

        self.read_to('eigvals(Ev)')
        orb_count = 0
        aos_read = False
        orb_start = False
        orb_done = False

        # Enter main loop to get molecular orbitals
        while orb_done == False and len(self.line) > 0:
            line = self.line.split()
            orbline = len(line) - 1
            orb_count += orbline
            #print line, orb_count, orbline

            if mo_min <= orb_count:
                # Read the relevant orbital info

                # Set up indices for loops later on
                if orb_start == False:
                    orb_index_a = mo_min - orb_count - 1
                else:
                    orb_index_a = -orbline
                if mo_max <= orb_count:
                    orb_index_b = mo_max - orb_count
                else:
                    orb_index_b = 0
       
                #print orb_index_a, orb_index_b
 
                # Add each MO to the list with an energy
                for i in range(1,len(line)):
                   x = MolOrb(float(line[i]))
                   self.mos.append(x)

                self.read_for(4) 
                        
                # Loop over each AO in this set of MOs
                for i in range(0,self.numorb):
                    line = self.line.split()
                    if len(line) == 0:
                        self.line = self.file.readline()
                        line = self.line.split()
    
                    # First orbital set only - store atomic orbital types
                    if aos_read == False:
                        at = int(self.line[5:9]) - 1
                        aotype = self.line[11:14].strip()
                        if 'xy2' in aotype: aotype = 'Dxy2'
                        self.atoms[at].add_aotype(aotype)

                        # Edit - store full AO name as aotype
                        #if 'P' in aotype: aotype = 'P'
                        #if 'D' in aotype or 'xy2' in aotype: aotype = 'D' 
                        
                        # Categorize all AOs by types from type file
                        a = []
                        for j in types:
                            #for listelem in j:
                            #   print listelem, aotype, listelem in aotype
                            if 'Atom' in j:
                                a.append(True) if (at >= (int(j[1])-1) and at < int(j[2])) else a.append(False)
                            elif 'Attype' in j:
                                a.append(True) if self.atoms[at].elem in j else a.append(False)
                            elif 'Orb' in j:
                                a.append(True) if any(listelem in aotype for listelem in j) else a.append(False)
                            elif 'Orbtype' in j:
                                a.append(True) if (any(listelem in aotype for listelem in j) and self.atoms[at].elem in j) else a.append(False)
                            elif 'Orbat' in j: 
                                a.append(True) if (at >= (int(j[1])-1) and at < int(j[2]) and any(listelem in aotype for listelem in j)) else a.append(False)
                            else:
                                 a.append(True)

                        self.atoms[at].add_ao(a)
                        #self.atoms[at].add_aotype(aotype) 
                        #print at, aotype, a

                        #if 'S' or 'P' or 'D' in type:
                        #    self.atoms[at].add_ao(type)
                        #else:
                        #    self.atoms[at].add_ao('Dxy2')

                    # Store the MO coefficients
                    # WARNING - ordering of i/j is reversed from previous scripts!!! ############
                    for j in range(orb_index_a,orb_index_b):
                        self.mos[j].add_coeff(line[j])

                    self.read_for(1)

                aos_read = True
                self.read_to('eigvals(Ev)','----')
                orb_start = True
                if mo_max <= orb_count:
                    orb_done = True

            else:
                # Skip to the next section
                self.read_to('eigvals(Ev)','----')

        # Compute MO character based on AO types
        for i in self.mos:
            for m in range(0,len(types)):
                i.char.append(0.0)
                ocycle = 0
                for j in self.atoms:
                    for k in j.ao:
                        if k[m]: i.char[m] += i.coeff[ocycle]**2
                        ocycle += 1
            #print i.char

        # Save that orbitals have been read
        if orb_count > 0:
            self.orbs_read = True
            #print "Orbitals read", orb_count
        else:
            print "Error: Orbitals not read"


    # Read electron configurations for CI
    def read_configs(self,types=[]):
        if self.states_read:
            self.file.seek(0)

        if not self.orbs_read:
            self.read_orbs(types)
        self.read_to('spin-adapted configurations')
        self.read_for(4)

        self.configs = []
        line = self.line.split()
        while len(self.line) > 5:
            ener = float(line[1])
            try:
                occ = int(self.line[53:57]) - 1
                vir = int(self.line[62:68]) - 1
            except ValueError:
                occ = 0
                vir = 0
            char = []
            achar = []
            bchar = []

            ochar = []
            vchar = []

            for i in range(0,len(self.mos[0].char)):
                char.append(self.mos[vir].char[i] - self.mos[occ].char[i])

                # achar = local character in orbitals included in type
                # bchar = local character in orbitals  outside of type
                achar.append(    min(self.mos[vir].char[i],self.mos[occ].char[i]))
                bchar.append(1 - max(self.mos[vir].char[i],self.mos[occ].char[i]))

                ochar.append(self.mos[occ].char[i])
                vchar.append(self.mos[vir].char[i])

                # RLG temporary edit for ligand clusters
                #char.append(self.mos[occ].char[i])
            self.configs.append(Config(occ,vir,char))
            self.configs[-1].achar = achar
            self.configs[-1].bchar = bchar
            self.configs[-1].ochar = ochar
            self.configs[-1].vchar = vchar
            self.configs[-1].change_energy(ener)

            #if len(self.configs) < 20:
            #    print char[0], achar[0], bchar[0], abs(char[0])+achar[0]+bchar[0]

            self.line = self.file.readline()
            line = self.line.split()

        # Save that configs have been read
        if len(self.configs) > 0:
            self.configs_read = True
            #print "CI configurations read"
        else:
            print "Error: CI configurations not read"


    # Read excited-state info
    def read_states(self,types=[]):
        # Make sure configs have been read
        if not self.configs_read and len(types) > 0:
            self.read_configs(types)

        self.read_to('CI trans.')
        self.read_for(3)

        self.states = [State(0.0,0.0)]

        # Read energies and oscillator strengths
        line = self.line.split()
        while len(line) > 2:
            num = int(line[0])
            try:
                e = float(line[1])
                o = float(line[4])
            except ValueError:
                e = float(line[2])
                o = float(line[5])
            self.states.append(State(e,o))
            if o > 1.e-5:
                pol = [float(self.line[50:60]),float(self.line[61:71]),float(self.line[72:82])]
            else:
                pol = [0.0,0.0,0.0]
            self.states[-1].calc_tdip(pol)
        
            self.line = self.file.readline()
            line = self.line.split()

        # Read CI coefficients
        if len(types) > 0 or self.configs_read == True:
            self.read_to('State')
            nstate = 0

            while 'State' in self.line:
                self.line = self.file.readline()

                ntypes = len(self.mos[0].char)
                char = []
                achar = []
                bchar = []

                # ochar & vchar for identifying character of occ/vir orbitals involved
                # Added for parametrization & auto-identification of states
                ochar = []
                vchar = []

                for i in range(0,ntypes): 
                    char.append(0.0)
                    achar.append(0.0)
                    bchar.append(0.0)

                    ochar.append(0.0)
                    vchar.append(0.0)

                while 'Config' in self.line:
                    line = self.line.split()
                    cnum = int(line[1])
                    coeff = float(line[3])
                    self.states[nstate].add_coeff(cnum, coeff)

                    if 'Config       1' in self.line:
                        while 'Config' in self.line:
                            self.line = self.file.readline()
                    else:
                        for i in range(0,ntypes):
                            char[i]  += coeff * self.configs[cnum-1].char[i]
                            achar[i] += coeff * self.configs[cnum-1].achar[i]
                            bchar[i] += coeff * self.configs[cnum-1].bchar[i]

                            ochar[i] += coeff * self.configs[cnum-1].ochar[i]
                            vchar[i] += coeff * self.configs[cnum-1].vchar[i]


                        self.line = self.file.readline()
                line = self.line.split()
                totcoeff = float(line[3])
                if totcoeff > 0:
                    for i in range(0,len(char)):
                        char[i] /= totcoeff
                        achar[i] /= totcoeff
                        bchar[i] /= totcoeff
                        ochar[i] /= totcoeff
                        vchar[i] /= totcoeff


                self.states[nstate].char = char
                self.states[nstate].achar = achar
                self.states[nstate].bchar = bchar
                self.states[nstate].ochar = ochar
                self.states[nstate].vchar = vchar

                #print nstate, char, ochar, vchar

                self.line = self.file.readline()
                self.line = self.file.readline()

                nstate += 1

            #for i in range(0,11):
            #    print self.states[i].energy, self.states[i].tdip

        # Save that states have been read
        if len(self.states) > 0:
            self.states_read = True
            #print "Excited states read"
        else:
            print "Error: Excited states not read"

    # Read Heat of Formation 
    def read_heat(self):
        if self.atoms_read or self.norbs_read or self.dipole_read or self.orbs_read or self.configs_read or self.states_read:
            self.file.seek(0)
 
        self.read_to('FINAL HEAT OF FORMATION')
        self.heat = float(self.line.split()[5]) if 'HEAT' in self.line else 0    
        
        # Save that heat of formation have been read
        if self.heat != 0: 
            self.heat_read = True
        else: 
            print "Error: Heat of formation not read"  



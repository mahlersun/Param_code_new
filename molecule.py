import sys
import string
import math
import constants as cnst
from properties import *

# Error function
def erf(x):
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

################################################################
# Parent molecule class                                        #
#                                                              #
# Use program-specific child classes (MopacMolecule, etc.) to  #
# input molecular information                                  #
#                                                              #
# Functions here handle post-processing of data and writing of #
# output files                                                 #
################################################################
class Molecule(object):
    def __init__(self, file):
        self.file = file
        self.line = self.file.readline()

        # Identify what has been read
        self.charge_read   = False
        self.atoms_read    = False
        self.norbs_read    = False
        self.dipole_read   = False
        self.orbs_read     = False
        self.configs_read  = False
        self.states_read   = False
        self.heat_read     = False # Newly added QS        
 
        # Dictionary of atomic symbols and numbers
        self.atomic_symbol = { 
            "H"  : 1 ,
            "He" : 2 ,
            "Li" : 3 ,
            "Be" : 4 ,
            "B"  : 5 ,
            "C"  : 6 ,
            "N"  : 7 ,
            "O"  : 8 ,
            "F"  : 9 ,
            "Ne" : 10,
            "Na" : 11,
            "Mg" : 12,
            "Al" : 13,
            "Si" : 14,
            "P"  : 15,
            "S"  : 16,
            "Cl" : 17,
            "Ar" : 18,
            "K"  : 19,
            "Ca" : 20,
            "Sc" : 21,
            "Ti" : 22,
            "V"  : 23,
            "Cr" : 24,
            "Mn" : 25,
            "Fe" : 26,
            "Co" : 27,
            "Ni" : 28,
            "Cu" : 29,
            "Zn" : 30,
            "Ga" : 31,
            "Ge" : 32,
            "As" : 33,
            "Se" : 34,
            "Br" : 35,
            "Kr" : 36,
            "Rb" : 37,
            "Sr" : 38,
            "Y"  : 39,
            "Zr" : 40,
            "Nb" : 41,
            "Mo" : 42,
            "Tc" : 43,
            "Ru" : 44,
            "Rh" : 45,
            "Pd" : 46,
            "Ag" : 47,
            "Cd" : 48,
            "In" : 49,
            "Sn" : 50,
            "Sb" : 51,
            "Te" : 52,
            "I"  : 53,
            "Xe" : 54,
            "Cs" : 55,
            "Ba" : 56,
            "La" : 57,
            "Ce" : 58,
            "Pr" : 59,
            "Nd" : 60,
            "Pm" : 61,
            "Sm" : 62,
            "Eu" : 63,
            "Gd" : 64,
            "Tb" : 65,
            "Dy" : 66,
            "Ho" : 67,
            "Er" : 68,
            "Tm" : 69,
            "Yb" : 70,
            "Lu" : 71,
            "Hf" : 72,
            "Ta" : 73,
            "W"  : 74,
            "Re" : 75,
            "Os" : 76,
            "Ir" : 77,
            "Pt" : 78,
            "Au" : 79,
            "Hg" : 80,
            "Tl" : 81,
            "Pb" : 82,
            "Bi" : 83,
            "Po" : 84,
            "At" : 85,
            "Rn" : 86
        }

    ################################################
    #                                              #
    #  General read functions                      #
    #                                              #
    ################################################

    # Read to a certain message in the file
    # Use a gibberish message2 to make it useless unless specified
    def read_to(self, message, message2='qwertyuiopasdf'):
        while message not in self.line and message2 not in self.line and len(self.line) > 0:
            self.line = self.file.readline()

    # Read a specific number of lines in the file
    def read_for(self,lines):
        for i in range(0,lines):
            self.line = self.file.readline()

    ################################################
    #                                              #
    #  Writing simple output files                 #
    #                                              #
    ################################################

    # Center a molecule relative to the origin
    def recenter(self):
        if not self.atoms_read:
            self.read_atoms(types)

        # Define a new origin
        max = [0.0,0.0,0.0]
        min = [0.0,0.0,0.0]
        origin = [0.0,0.0,0.0]
        for i in range(0,3):
            max[i] = max(self.atoms, key = lambda a: a.coord[i])
            min[i] = min(self.atoms, key = lambda a: a.coord[i])
            origin[i] = (max[i]+min[i])/2

        print origin
        # Displace atoms
        for a in self.atoms:
            for i in range(0,3):
                a.coord[i] -= origin[i]

        '''
        max = [xyz[0][1],xyz[0][2],xyz[0][3]]
        min = [xyz[0][1],xyz[0][2],xyz[0][3]]
        for i in range(1,numat):
            for j in range(0,3):
                if xyz[i][j+1] > max[j]:
                    max[j] = xyz[i][j+1]
                elif xyz[i][j+1] < min[j]:
                    min[j] = xyz[i][j+1]
        origin = [(max[0]+min[0])/2,(max[1]+min[1])/2,(max[2]+min[2])/2]
        
        for i in range(1,numat):
            for j in range(0,3):
                xyz[i][j+1] -= origin[j]
        '''

    # Using a template, write an input file using coordinates read in
    def write_input(self, outfilename, template):
        if not self.charge_read:
            self.read_charge()
        if not self.atoms_read:
            self.read_atoms()

        # A few default options for common outputs
        if 'def' in template:
            if 'xyz' in template:
                out = open(outfilename+'.xyz','w')
                out.write(str(len(self.atoms)) + '\n\n')
                for a in self.atoms:
                    out.write(string.rjust(a.elem,4))
                    for j in range(0,3):
                        out.write(string.rjust('%.5f'%a.coord[j],12))
                    out.write('\n')
                out.write('\n')
                out.close()
            elif 'mopac' in template:
                out = open(outfilename+'.dat','w')
                out.write('INDO CHARGE='+str(self.charge)+' RCI\nTitle\nINDO\n')
                for a in self.atoms:
                    out.write(string.rjust(a.elem,4))
                    for j in range(0,3):
                        out.write(string.rjust('%.5f'%a.coord[j],12) + ' 0')
                    out.write('\n')
                out.write('\n')
                out.close()
        # Generic output using templates
        # Template may include NAT, CHRG, TYPE, ELEM, XYZ, XXX, YYY, ZZZ as needed
        # Append suffix of template onto outfilename
        else:
            temp = open(template,'r')
            out = open(outfilename + template[template.rfind('.'):],'w')

            tline = temp.readline()
            while len(tline) > 0:
                # Single replacements
                if 'NAT' in tline:
                    tline = tline[:tline.find('NAT')] + str(len(self.atoms)) + tline[tline.find('NAT')+3:]
                if 'CHRG' in tline:
                    tline = tline[:tline.find('CHRG')] + str(self.charge) + tline[tline.find('CHRG')+4:]
                if 'TITLE' in tline:
                    tline = tline[:tline.find('TITLE')] + outfilename + tline[tline.find('TITLE')+5:]

                # Use line as a template for printing all atoms
                if 'ELEM' in tline or 'XXX' in tline:
                    for a in self.atoms:
                        line = tline
                        if 'ELEM' in tline:
                            line = line[:line.find('ELEM')] + string.rjust(a.elem,4) + line[line.find('ELEM')+4:]
                        if 'XYZ' in tline:
                            line = line[:line.find('XYZ')] + string.rjust('%.5f'%a.coord[0],12) + string.rjust('%.5f'%a.coord[1],12) + string.rjust('%.5f'%a.coord[2],12) + line[line.find('XYZ')+3:]
                        if 'XXX' in tline:
                            line = line[:line.find('XXX')] + string.rjust('%.5f'%a.coord[0],12) + line[line.find('XXX')+3:]
                        if 'YYY' in tline:
                            line = line[:line.find('YYY')] + string.rjust('%.5f'%a.coord[1],12) + line[line.find('YYY')+3:]
                        if 'ZZZ' in tline:
                            line = line[:line.find('ZZZ')] + string.rjust('%.5f'%a.coord[2],12) + line[line.find('ZZZ')+3:]
                        out.write(line)
                # Loop over atom types (use for ADF basis sets)
                elif 'TYPE' in tline:
                    for t in self.at_types:
                        line = tline
                        while 'TYPE' in tline:
                            line = line[:line.find('TYPE')] + string.rjust(t,4) + line[line.find('TYPE')+4:]
                        out.write(line)
                # Print line with edits as needed
                else:
                    out.write(tline)
                tline = temp.readline()
            out.close()

    # Write vibrational modes and other info needed for Raman calculations
    def write_vibs(self, outfilename, minfreq=500, maxfreq=2000):
        return

    def write_orbs(self, outfilename, types=[]):
        if not self.orbs_read:
            self.read_orbs(types)

        orb = open(outfilename+'.orb','w')
        orb.write("  Num    Energy")

        for j in types:
            orb.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        orb.write('\n')

        for i in range(0,len(self.mos)):
           orb.write(string.rjust(str(i+1), 5) + string.rjust("%.4f" % self.mos[i].energy, 10))
           for j in range(0,len(types)):
               orb.write(string.rjust("%.4f" % self.mos[i].char[j], 12))
           orb.write('\n')
        orb.close()

        #obin = open(outfilename+'.orb_bin','w')
        # Finish here??

    def write_omo(self, outfilename):
        if not self.orbs_read:
            self.read_orbs([])
        omo = open(outfilename+'.omo','w')
        omo.write(string.rjust(str(self.numorb),6) + '\n')

        nstart = 0
        perline = 6
        
        while nstart < self.numorb:
            nend = nstart + perline
            if nend > self.numorb:
                nend = self.numorb
            for i in range(0,self.numorb):
                for j in range(nstart,nend):
                    omo.write(string.rjust('%.6f'%(self.mos[j].coeff[i]),10))
                omo.write('\n')
            nstart += perline

        omo.close()

    # Write oscillator strength info
    def write_osc(self, outfilename, types=[]):
        if not self.states_read:
            self.read_states(types)
        
        osc = open(outfilename+'.osc','w')
        osc.write("  Num    Energy     Osc Str")

        for j in types:
            osc.write(string.rjust(j[0][:3] + ' ' + j[1], 12))

        osc.write('\n')

        nstate = len(self.states)

        for j in range(1,nstate):
            osc.write(string.rjust(str(j+1),5) + string.rjust("%.4f" % self.states[j].energy,10) + string.rjust("%.5f" % self.states[j].osc,12))

            for i in range(0,len(types)):
                osc.write(string.rjust("%.4f" % self.states[j].char[i], 12))

            osc.write('\n')

        osc.close()

    # Write excited-state analysis
    def write_exc(self, outfilename, types=[]):
        if not self.configs_read:
            self.read_mos(types)
            self.read_configs(types)
            self.read_states(types)
        if not self.states_read:
            self.read_states(types)

        for state in self.states:
            state.find_energy_shift(self.configs, self.mos)
            state.find_collectivity()
        self.find_superatom_char(types)

        nstate = len(self.states)

        try: 
            dip = open(outfilename+'.dip','r')
            dline = dip.readline()
            for j in range(1,nstate):
                dline = dip.readline()
                line = dline.split()
                self.states[j].dipadd = float(line[3])/100.
            dipadd = True
            dip.close()
        except IOError:
            dipadd = False

        exc = open(outfilename+'.exc','w')
        ex2 = open(outfilename+'.ex2','w')
        exc.write("#Num    Energy   Osc Str Superatom   Collect CoupRange Interband")
        ex2.write("#Num    Energy   Osc Str Superatom   Collect CoupRange Interband")
        if dipadd == True:
            exc.write("   Dip Add")
            ex2.write("   Dip Add")
        exc.write('\n')
        ex2.write('\n')

        for j in range(1,nstate):
            exc.write(string.rjust(str(j+1),4) + string.rjust("%.4f" % self.states[j].energy,10) + string.rjust("%.5f" % self.states[j].osc,10))
            exc.write(string.rjust("%.4f" % self.states[j].superatom_char,10) + string.rjust("%.4f" % self.states[j].collectivity,10)  )
            exc.write(string.rjust("%.4f" % self.states[j].weighted_std,10)   + string.rjust("%.4f" % self.states[j].interband_char,10))
            if dipadd == True:
                exc.write(string.rjust("%.4f" % self.states[j].dipadd,10))
            exc.write('\n')
        exc.write('\n\n\n')


#        maxnume = min(200,len(self.states))-1
#        e_max = min(self.states[maxnume].energy,7.0)

#        self.states.sort(key=lambda state: state.weighted_std, reverse=True)
        # Print only states with non-negligible transition dipoles

#        ex2.write("#Num    Energy   Osc Str Superatom   Collect CoupRange Interband")
#        if dipadd == True:
#            ex2.write("   Dip Add")
#        ex2.write('\n')
        for j in range(1,nstate):
            if self.states[j].osc > 1e-2:# and self.states[j].superatom_char > 0.1 and self.states[j].collectivity > 2:
                ex2.write(string.rjust(str(j+1),4) + string.rjust("%.4f" % self.states[j].energy,10) + string.rjust("%.5f" % self.states[j].osc,10))
                ex2.write(string.rjust("%.4f" % self.states[j].superatom_char,10) + string.rjust("%.4f" % self.states[j].collectivity,10)  )
                ex2.write(string.rjust("%.4f" % self.states[j].weighted_std,10)   + string.rjust("%.4f" % self.states[j].interband_char,10))
                if dipadd == True:
                    ex2.write(string.rjust("%.4f" % self.states[j].dipadd,10))
                ex2.write('\n')

        exc.close()





    # Compute the absorption spectrum
    def write_sigma(self, outfilename, types=[], max_e=8.0, e_step=0.02, gamma=0.2):
        # Make sure states have been read
        if not self.states_read:
            self.read_states(types)

        step_count = int(max_e/e_step) + 1
        sigma = []
        for i in range(0,step_count):
            sigma.append([i*e_step, 0.0])
            for j in range(0,len(types)*2):
                sigma[i].append(0.0)

        for j in self.states:
            # Add the cross-section to the appropriate bins
            if j.osc > 1e-6:
                for i in range(0,step_count):
                    current_e = sigma[i][0]
        
                    # Use Lorentzian line shape
                    phi = gamma/(((current_e - j.energy)**2 + gamma**2)*math.pi)
        
                    sigma[i][1] += phi * j.osc

                    for k in range(0,len(types)):
                        if j.char[k] > 0.0: 
                            sigma[i][2 + k*2] += phi * j.osc * j.char[k]
                        else:
                            sigma[i][3 + k*2] -= phi * j.osc * j.char[k]

        sig = open(outfilename+'.sigma','w')

        sig.write('  Energy        Tot Abs')
        for j in types:
            sig.write(string.rjust('To ' + j[0][:3] + ' ' + j[1], 15) + string.rjust('From ' + j[0][:3] + ' ' + j[1], 15))
        sig.write('\n')

        for i in range(0,step_count):
            sig.write(string.rjust("%.4f"%sigma[i][0],8)) 
            for j in sigma[i][1:]:
                sig.write(string.rjust("%.6f" % j ,15))
            sig.write('\n')
        sig.close()

    def compute_alpha(self, omega=0.0, gamma=0.1088j, states=0, pr=False, gamma2=0.0, types=[], usefreq=False, modefreq=0.0):
        # Make sure states have been read
        #if not self.orbs_read:
        #    self.read_orbs([])
        if not self.states_read:
            self.read_states(types)

        alpha = [[[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]]
        #for i in range(0,len(types)*3):
        #    alpha.append([[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],
        #                  [complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]])

        if states == 0 or states >= len(self.states):
            states = len(self.states) - 1

        g_state = gamma
        e_diff = modefreq / cnst.ev2cm

        for k in range(1,states+1):
            s = self.states[k]
            #print s.energy, e_diff, omega, gamma, s.energy + omega + g_state, s.energy - e_diff - omega - g_state, s.energy - omega - g_state, s.energy - e_diff + omega + g_state
            #print k, s.char[0], s.achar[0],s.bchar[0],s.char[0]+s.achar[0]+s.bchar[0]

            # Compute state contribution to alpha
            #a_state = [[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)],[complex(0.0,0.0),complex(0.0,0.0),complex(0.0,0.0)]]
            
            # Compute a state-specific gamma based on two input gamma values and excitation character
            if gamma2 != 0.0:
                if abs(s.char[0]) > 0.9999:
                    g_state = gamma2
                elif abs(s.char[0]) < 0.0001:
                    g_state = gamma
                else:
                    # Linear scaling - huge change in gamma for states with predominantly one character
                    #g_state = abs(s.char[0])*gamma2 + (1-abs(s.char[0]))*gamma

                    # Quadratic scaling
                    #g_state = gamma + abs(s.char[0])**2 * (gamma2 - gamma)

                    # Gaussian scaling
                    gscl = 4
                    g_state = (gamma + gamma2)/2 + erf(gscl*(abs(s.char[0])-0.5))/erf(gscl*0.5)*(gamma2 - gamma)/2


            # Following Reimers code, convert to cm-1 here
            if usefreq == False:
                a_denom_1 = (s.energy - g_state - omega) * cnst.ev2cm
                a_denom_2 = (s.energy + g_state + omega) * cnst.ev2cm
                #print '1',a_denom_1/cnst.ev2cm, a_denom_2/cnst.ev2cm
            else:
                # Explicitly include the difference in energies between the initial and final vibronic states
                e_diff = modefreq / cnst.ev2cm

                # Version 1: Standard Kramers-Heisenberg formula
                #a_denom_1 = (s.energy          - omega - g_state) * cnst.ev2cm
                #a_denom_2 = (s.energy - e_diff + omega + g_state) * cnst.ev2cm  #Stokes line where (s.energy - e_diff) = transition energy from s to final state

                #alpha[0][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph

                # Version 2: Standard Kramers-Heisenberg formula with states reversed
                a_denom_1 = (s.energy          + omega + g_state) * cnst.ev2cm
                a_denom_2 = (s.energy - e_diff - omega - g_state) * cnst.ev2cm  #Stokes line where (s.energy - e_diff) = transition energy from s to final state
                #print '2',a_denom_1/cnst.ev2cm, a_denom_2/cnst.ev2cm

            for i in range(0,3):
                for j in range(0,3):
                    a_num = s.tdip[i]*s.tdip[j]
                    alpha[0][i][j] += (a_num/a_denom_1 + a_num/a_denom_2) * cnst.alph



        if pr == True:
            for k in range(0,len(alpha)):
                print 'Total', string.rjust('%.4f'%abs((alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0),10), string.rjust('%.4f'%((alpha[k][0][0].real + alpha[k][1][1].real + alpha[k][2][2].real) / 3.0),10), string.rjust('%.4f'%((alpha[k][0][0].imag + alpha[k][1][1].imag + alpha[k][2][2].imag) / 3.0),10)

        # Orientationally average alpha
        #a_or = complex(0.0,0.0)
        #a_or +=  (alpha[0][0] + alpha[1][1] + alpha[2][2]) / 3.0
        #print a_or, abs(a_or)

        return alpha

    def write_alpha(self, outfilename, omega=0.0, gamma=0.1088j, states=0, gamma2=0.0, types=[]):
        alpha = self.compute_alpha(omega, gamma, states, True, gamma2, types)

        #a_or =  (alpha[0][0][0] + alpha[0][1][1] + alpha[0][2][2]) / 3.0

        al_out = open(outfilename + '.alpha', 'w')
        al_out.write('Omega = ' + '%.4f'%omega + ' eV   Gamma = ' + '%.4f'%gamma.imag + ' eV')
        if gamma2 != 0.0:
            al_out.write('Gamma2 = ' + '%.4f'%gamma2.imag + ' eV')
        al_out.write('\n\n')

        for k in range(0,len(alpha)):
            a_or =  (alpha[k][0][0] + alpha[k][1][1] + alpha[k][2][2]) / 3.0

            al_out.write('Real polarizability (au)\n')
            print 'Real polarizability (au)'
            for i in range(0,3):
                al_out.write(string.rjust('%.4f'%alpha[k][i][0].real,12) + string.rjust('%.4f'%alpha[k][i][1].real,12) + string.rjust('%.4f'%alpha[k][i][2].real,12) + '\n')
                print string.rjust('%.4f'%alpha[k][i][0].real,12) + string.rjust('%.4f'%alpha[k][i][1].real,12) + string.rjust('%.4f'%alpha[k][i][2].real,12)
            al_out.write('\n')
            al_out.write('Imag polarizability (au)\n')
            print 'Imag polarizability (au)'
            for i in range(0,3):
                al_out.write(string.rjust('%.4f'%alpha[k][i][0].imag,12) + string.rjust('%.4f'%alpha[k][i][1].imag,12) + string.rjust('%.4f'%alpha[k][i][2].imag,12) + '\n')
                print string.rjust('%.4f'%alpha[k][i][0].imag,12) + string.rjust('%.4f'%alpha[k][i][1].imag,12) + string.rjust('%.4f'%alpha[k][i][2].imag,12)
            al_out.write('\n')

        
            al_out.write('Orientationally averaged polarizability (au)\n')
            al_out.write(string.rjust('%.4f'%a_or.real,12) + string.rjust('%.4f'%a_or.imag,12) + string.rjust('%.4f'%abs(a_or),12) + '\n\n\n')

        al_out.close()

    # Find the superatomic character
    # Based on Adam's PEW script
    def find_superatom_char(self, types, cutoff=0.4):
        # Check if type information has been calculated
        t = False
        for i in range(0,len(types)):
            type = types[i]
            if 'SUPER' in type:
                t_ind = i
                t = True
        if t == False:
            types.append(['Orbtype','Ag','S','P'])
            t_ind = len(types) - 1
            self.read_orbs(types)
            self.read_configs(types)
            self.read_states(types)
        print 't_ind',t_ind

        # Determine whether each orbital is superatomic
        for mo in self.mos:
            if mo.char[t_ind] > cutoff:
                mo.super = True
            else:
                mo.super = False
            #print mo.super,mo.char[t_ind]

        # Determine the superatomic character of each excited state
        # Depends only on occupied orbital characters
        for exc in self.states:
            exc.superatom_char = 0.0
            exc.interband_char = 0.0
            for conf in exc.coeff:
                # conf[0] = configuration number; conf[1] = weight
                for occ in self.configs[conf[0]-1].occ:
                    exc.interband_char += conf[1]*self.mos[occ].char[t_ind]
                    if self.mos[occ].super == True:
                        exc.superatom_char += conf[1]





import math

fosc_fact = 0.0038
#alph = 0.000029376024  # Used in polarizability calculation

# Energy units
ev2cm  = 8065.73
au2ev  = 27.2107
au2cm  = 219474.63
au2ang = 0.5291772
debye  = 4.80294

bohr2ang = 0.5291772

# Polarizability unit conversion
alph = au2cm/au2ang**2/debye**2

planck = 6.6260755E-34
boltz = 1.3806580E-23
speed = 2.99792458E+8
epsilon0 = 8.8541878E-12
avogadro = 6.02214199E+23
amu = 1.0/avogadro*1E-03
temper = 300.0
exparg = planck*speed*100.0/boltz
conver = 2.0*planck*math.pi**2/speed*1.0E-40/amu
width = 20.0
scale = 1E+34

afac = width/(2.0*math.pi)
bfac = width/2.0

# Atomic masses
mass_dict = {
    'H'  :   1.008, 
    'He' :   4.0026,
    'Li' :   6.94,
    'Be' :   9.0122,
    'B'  :  10.81,
    'C'  :  12.011,
    'N'  :  14.007,
    'O'  :  15.999,
    'F'  :  18.998,
    'Na' :  22.989,
    'Mg' :  24.305,
    'Al' :  26.981,
    'Si' :  28.085,
    'P'  :  30.973,
    'S'  :  32.06,
    'Cl' :  35.45,
    'K'  :  39.0983,
    'Ca' :  40.078,
    'Fe' :  55.845,
    'Co' :  58.933,
    'Ni' :  58.6934,
    'Cu' :  63.546,
    'Zn' :  65.38,
    'As' :  74.9216,
    'Se' :  78.96,
    'Br' :  79.904,
    'Ru' : 101.07,
    'Rh' : 102.90,
    'Pd' : 106.42,
    'Ag' : 107.8682,
    'Cd' : 112.411,
    'In' : 114.818,
    'Sn' : 118.71,
    'I'  : 126.90,
    'Ir' : 192.217,
    'Pt' : 195.084,
    'Au' : 196.96,
    'Hg' : 200.59
}






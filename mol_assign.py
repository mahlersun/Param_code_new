from mopac_mol import *
from adf_mol   import *
from qchem_mol import *
from man_mol   import *
from gauss_mol import *

def open_mol(filename, prog = ' '):
    # Open input and output files
    try:
        file = open(filename,'r')
    except IOError:
        try:
            file = open(filename+'.out','r')
        except IOError:
            file = open(filename+'.log','r')   
 
    # Check what program this input is associated with and open as the appropriate class
    if prog == ' ':
        fline = file.read(5000).lower()
        #print fline
        file.seek(0)
        if 'mopac: public domain version' in fline or 'mopac2016' in fline:
            prog = 'mopac'
        elif 'amsterdam density functional' in fline:
            prog = 'adf'
        elif 'q-chem, inc., pleasanton, ca' in fline:
            prog = 'qchem'
        elif 'entering gaussian' in fline:
            prog = 'gauss'
        else:
            # Assume manual input
            prog = 'man'
    
    if prog == 'mopac':
        out = MopacMolecule(file)
    elif prog == 'adf':
        out =   AdfMolecule(file)
    elif prog == 'qchem':
        out =  QchemMolecule(file)
    elif prog == 'gauss':
        out =  GaussMolecule(file)
    elif prog == 'man':
        out =   ManMolecule(file)
    else:
        sys.exit("Program"+prog+"not recognized")

    return out, prog

def read_types(tfilename):
    types = []
    try:
        tfile = open(tfilename,'r')
        for line in tfile.readlines():
            print 'Type   %s' % line,
            types.append(line.split())
    except IOError:
        print '%s not found; no types assigned' % tfilename

    return types



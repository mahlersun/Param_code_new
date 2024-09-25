#!/home/gieseking/bin/miniconda2/envs/python2sci/bin/python

#-------------------------------------------------------------------------
#
# For a list of MOPAC PM7 file names:
#   1. Call MOPAC to compute 
#   2. Read in reference data from comment line of .dat file
#   3. Compute and print error statistics for all files
#      
# Last edited by: Qiwei Sun 12/2021
#-------------------------------------------------------------------------

import string,sys,math
import constants as cnst
from mol_assign  import *
from ref_data    import *
from param  import *
import getopt

tfilename = 'types'
types = []
prog = ' '
processors = 1

helpfile = """
param_pm.py -i <inputfile> -o <outputfile>
    -i   Give the name of a .inp file containing the parametrization procedure. The format of the file is:

         ParamList      # List atom number, then parameter name; follow by allowed ranges.
           47 betas   -1.5   -4.0
           47 betad  -25.0  -50.0
         End
         
         ParamVal       # List one parameter set per line, in same order as ParamList
           -1.5  -15.0
         End
         
         InitParamFile  # INDO parameter file to edit
           starting_point.txt
         End
         
         MolList        # List all molecules, one per line
           Ag2_0_01
           Ag3_1_01
         End
         
         Procedure      # Parametrization procedure 
           ReadError       True       # Read previously computed errors?
             ErrorFile     param.log  # File to read/write errors
             ClearOld      False      # Overwrite old data or append new to existing file?
           GenParam        False      # Generate new parameters randomly?
           GenParamSparse  True       # Generate new parameters randomly, enforcing minimum
                                      # spacing between parameter sets?
             RandCnt       50         # Number of random parameter sets to generate
             SparseDist    0.15       # Minimum spacing of sparse parameter sets
           NeighborIter    4          # Iterations of finding neighbors and generating new
                                      # nearby random parameter sets
         End
"""

# Parse input options
try:
    options, remainder = getopt.getopt(sys.argv[1:],"hi:o:t:p:",['--help','--input=','--output=','--types=','--ncpus='])
except getopt.GetoptError as err:
    print str(err)
    print helpfile
    sys.exit(2)

if remainder != []:
    print "Error: options not read - %r" % remainder
    print helpfile
    sys.exit(2)


for opt, arg in options:
    if opt in ('-h','--help'):
        print helpfile
        sys.exit(2)
    elif opt in ('-i','--input'):
        #if '.' in arg:
        #    arg = arg[:arg.rfind('.')]
        infilename = arg
        outfilename = arg
    elif opt in ('-o','--output'):
        if '.' in arg:
            arg = arg[:arg.rfind('.')]
        outfilename = arg
    elif opt in ('-t','--type'):
        tfilename = arg
    elif opt in ('-p','--ncpus'):
        processors = int(arg)         

# Set up defaults for procedural options
param_list       = []
param_range      = []
param_val_list   = []
init_param_file  = 'startpoint.txt' # change
mol_list         = []
read_error       = False
error_file       = 'param.log'
clear_old        = False
gen_param        = False
gen_param_sparse = False
rand_count       = 5
sparse_dist      = 0.15
neighbor_iter    = 0
graph_error      = False

# Read procedural options from input file
inp = open(infilename,'r')
iline = inp.readline()
while len(iline) > 0:
    if 'paramlist' in iline.lower():
        iline = inp.readline()
        while 'end' not in iline.lower() and len(iline) > 0:
            line = iline.split()
            param_list.append([int(line[0]),line[1]])
            if len(line) > 3:
                param_range.append([float(line[2]),float(line[3])])
            else:
                param_range.append([])
            iline = inp.readline()
        if 'end' not in iline.lower():
            print "Error: No end to ParamList input block"
        #print 'param_list', param_list, param_range
    if 'paramval' in iline.lower():
        iline = inp.readline()
        while 'end' not in iline.lower() and len(iline) > 0:
            line = iline.split()
            param_val = []
            for par in line:
                param_val.append(float(par))
            param_val_list.append(param_val)
            iline = inp.readline()
        if 'end' not in iline.lower():
            print "Error: No end to ParamVal input block"
        #print 'param_val', param_val_list
    if 'mollist' in iline.lower():
        iline = inp.readline()
        while 'end' not in iline.lower() and len(iline) > 0:
            mol_list.append(iline.split()[0])
            iline = inp.readline()
        if 'end' not in iline.lower():
            print "Error: No end to MolList input block"
        #print 'mol_list', mol_list
    if 'procedure' in iline.lower():
        iline = inp.readline()
        while 'end' not in iline.lower() and len(iline) > 0:
            line = iline.split()
            if   line[0].lower() == 'initparamfile':
                init_param_file   = line[1]
            elif line[0].lower() == 'readerror':
                if line[1].lower() in ['true','t','yes','y','1']: read_error        = True
            elif line[0].lower() == 'errorfile':
                error_file        = line[1]
            elif line[0].lower() == 'clearold': 
                if line[1].lower() in ['true','t','yes','y','1']: clear_old         = True
            elif line[0].lower() == 'genparam':
                if line[1].lower() in ['true','t','yes','y','1']: gen_param         = True
            elif line[0].lower() == 'genparamsparse':
                if line[1].lower() in ['true','t','yes','y','1']: gen_param_sparse  = True
            elif line[0].lower() == 'randcount':
                rand_count        = int(line[1])
            elif line[0].lower() == 'sparsedist':
                sparse_dist       = float(line[1])
            elif line[0].lower() == 'neighboriter':
                neighbor_iter     = int(line[1])
            elif line[0].lower() == 'grapherror':
                if line[1].lower() in ['true','t','yes','y','1']: graph_error       = True
            elif len(line) > 0:
                print "Error: Option not recognized  ", line
            iline = inp.readline()
        if 'end' not in iline.lower():
            print "Error: No end to Procedure input block"
        #print 'procedure', init_param_file, read_error, error_file, clear_old, gen_param, gen_param_sparse, rand_count, sparse_dist, neighbor_iter, graph_error
    iline = inp.readline()

if len(param_list) == 0:
    print "Error: No tuable parameters read"

# Algorithm Starts Here
par = Param_Set(mol_list, param_list, init_param_file, param_val_list, read_error=read_error, param_range=param_range, clear_old=clear_old, logfile=error_file, processors=processors)

if gen_param:
    par.gen_params(rand_count, param_range)

if gen_param_sparse:
    par.gen_params_sparse(rand_count, sparse_dist, param_range)

for k in range(0,neighbor_iter):
    count = len(par.params) # number of paramters sets
    near = max(3, int(math.sqrt(count))) # the greater one between 3 and (sqrt of count)
    new = min(len(param_list)*5, int(math.sqrt(count))) # least one between the dimensions of parameter vector and the sqrt of count
   # new = min(len(param_list)*5, int(math.sqrt(count))/2)
    loc_min_list, distances = par.find_local_min(near)
    d_last = new if new < len(distances[0]) else -1
    for i in loc_min_list:
        par.gen_params_near(new, i, distances[i][d_last], par.param_range)
count = len(par.params)
near = max(3, int(math.sqrt(count)))
#print 'count', count, near
loc_min_list, distances = par.find_local_min(near)

print "Total count: ", len(distances)

if graph_error:
    par.plot_error_step()



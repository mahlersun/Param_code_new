###############################################################
# Class to generate and modify parameters for                 #
# reparametrization.                                          #
#                                                             #
# Class Param_Set:                                            #
#   Contains data for a series of parameter sets              #
#   __init__: Set up initial values                           #
#                                                             #
#   gen_params: Given a list of parameters (param_list) and   #
#       allowed ranges (param_ranges), use random numbers to  #
#       generate parameter values to store in param_val_list  #
#                                                             #
#                                                             #
# Class Param:                                                #
#   Contains data for one parameter set                       #
#   __init__: Set up initial values                           #
#                                                             #
#   read_orbs: For indo_dens; read basis set                  #
#                                                             #
#   edit_param: Given a list of parameters (param_list) and   #
#       values(param_val), edit INDO1S.par to contain         #
#       the new values                                        #
#                                                             #
#   run_mopac: Run MOPAC for all input files in mol_list      #
#                                                             #
#   param_error: For a given set of parameters, compute and   #
#        store the RMS error for the given set of reference   #
#        data                                                 #
#                                                             #
###############################################################

import string
import os
import random
from ref_data    import *
from mol_assign  import *
from relative_energy import *
import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
import time
#%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d    import Axes3D
from sklearn.preprocessing   import PolynomialFeatures
from sklearn.linear_model    import LinearRegression
from sklearn.metrics         import mean_squared_error, r2_score
#from sklearn.svm             import SVR
from sklearn.kernel_ridge    import KernelRidge
from sklearn.neighbors       import KNeighborsRegressor, NearestNeighbors

class Param_Set(object):
    def __init__(self, mol_list, param_list, init_param_file, param_val_list=[],logfile='param.log', read_error=False, param_range=[], clear_old=False, processors=1):
        self.mol_list = mol_list
        if read_error and clear_old:
            print "Error: Cannot read and clear errors. Defaulting to read without clearing"
            self.log = open(logfile,'a+')
        elif read_error and not clear_old:
            self.log = open(logfile,'a+')
        elif clear_old:
            self.log = open(logfile,'w')
        else:
             self.log = open(logfile,'a')
        self.ref_data = Ref_Set(mol_list)
        self.ref_rele = Ref_Rel_E(mol_list, self.ref_data)
        self.param_range = param_range
        self.params_rescale = []
        self.processors = processors

        # Parameter file to be edited (name cannot be INDO1S.par)
        self.init_param_file = init_param_file

        # List of parameters to tune, sorted by atomic number and parameter order
        self.param_list = param_list

        # List of parameter sets
        self.params = []
        if read_error == True:
            self.read_error()
        for param_val in param_val_list:
            self.params.append(Param(self.param_list, self.init_param_file, param_val))
            self.params[-1].param_error(self.mol_list, self.ref_data, self.log, self.processors, self.ref_rele)


    # Rescale parameters to intervals of [0,1] based on allowed ranges
    def param_rescale(self):
        self.params_rescale = []
        for par in self.params:
            self.params_rescale.append([])
            for i in range(0,len(self.param_range)):
                self.params_rescale[-1].append((par.param_val[i] - self.param_range[i][0])/(self.param_range[i][1] - self.param_range[i][0]))
            #print par.param_val, self.params_rescale[-1] 


    # Find nearest neighbors for parameters
    def find_neighbors(self, nneigh):
        self.param_rescale()
        nbrs = NearestNeighbors(n_neighbors=nneigh).fit(self.params_rescale)
        distances, indices = nbrs.kneighbors(self.params_rescale)
        #print distances 
        return distances, indices


    # Find parameter sets with lower error than nearest neighbors
    def find_local_min(self, nneigh):
        loc_min_list = []
        distances, indices = self.find_neighbors(nneigh)
        err_list = []
        for i in range(0,len(self.params)):
            err_list.append(self.params[i].rms)
        cutoff = 3 * min(err_list) # Use an arbitrary y = n for cut off, y may vary later
        for i in range(0,len(self.params)):
            loc_min = True 
            err_curr = self.params[i].rms
            for j in indices[i][1:]:
                if err_curr > self.params[j].rms:
                    loc_min = False
            if loc_min == True:
                if err_curr < cutoff:  
                    loc_min_list.append(i) 

        # Introduce cross-breeding        
        param_cross = []
        for i in range(1,len(loc_min_list)):
            f = loc_min_list[i]
            m = loc_min_list[i-1]
            param_temp = []
            for j in range(0, len(self.param_list)):
                param_temp.append((self.params[f].param_val[j] + self.params[m].param_val[j]) / 2)
            param_cross.append(param_temp)
        for i in range(0, len(param_cross)):
            self.params.append(Param(self.param_list, self.init_param_file, param_cross[i]))
            self.params[-1].param_error(self.mol_list, self.ref_data, self.log, self.processors, self.ref_rele)
            print 'RMS of the hybrid is ', '%.4f'%self.params[-1].rms 
           # loc_min_list.append((len(self.params) - 1))
      # print loc_min_list
        print 'Local minima'
        for i in loc_min_list:
            print i, '%.4f'%self.params[i].rms, '%.4f'%distances[i][-1]
        #    self.gen_params_near(5, i, distances[i][-1], self.param_range)
        return loc_min_list, distances


    # Read previously computed errors from a file
    # WARNING: No error checking of molecules or parameter list 
    def read_error(self):
        len_init = len(self.params)
        lline = self.log.readline()
        while 'Param' in lline:
            lline = self.log.readline()
            param_val = self.log.readline().split()
            for i in range(0,len(param_val)):
                param_val[i] = float(param_val[i])
            self.params.append(Param(self.param_list, self.init_param_file, param_val))
            while 'Weighted' not in lline and len(lline) > 0:
                lline = self.log.readline()
            line = lline.split()
            self.params[-1].mse = float(line[-1])
            lline = self.log.readline()
            line = lline.split()
            self.params[-1].mae = float(line[-1])
            lline = self.log.readline()
            line = lline.split()
            self.params[-1].rms = float(line[-1])
            lline = self.log.readline()
            lline = self.log.readline()
            #print 'param_val', param_val
            #print 'Error', self.params[-1].mse, self.params[-1].mae, self.params[-1].rms
        print (len(self.params) - len_init), 'parameter sets read from file'

    # Generate random parameters based on a list
    def gen_params(self, number, param_ranges):
        print 'Generating ', number, ' random parameter sets'
        for i in range(0,number):
            param_val = []
            for j in range(0, len(self.param_list)):
                a = random.random()
                param_val.append(param_ranges[j][0] + a * (param_ranges[j][1] - param_ranges[j][0]))
                #print a, param[-1]
            self.params.append(Param(self.param_list, self.init_param_file, param_val))
            #print i, self.param_val_list
            self.params[-1].param_error(self.mol_list, self.ref_data, self.log, self.processors, self.ref_rele)


    # Generate parameters far from existing points
    def gen_params_sparse(self, number, distmin, param_ranges):
        print 'Generating ', number, ' parameter sets ', '%.4f'%distmin, ' from existing sets'
        if len(self.params_rescale) < len(self.params):
            self.param_rescale()

        params_added, params_tried = 0, 0
        while params_added < number and params_tried < number*100:
            param_val_rescale = []
            valid = True
            for j in range(0, len(self.param_list)):
                param_val_rescale.append(random.random())
            # Check if far enough from other parameter sets
            for i in range(0,len(self.params_rescale)):
                dist = 0.0
                for j in range(0, len(self.param_list)):
                    dist += (param_val_rescale[j] - self.params_rescale[i][j])**2
                dist = math.sqrt(dist)
                if dist < distmin:
                    valid = False
            if valid == True:
                param_val = []
                for j in range(0, len(self.param_list)):
                    param_val.append(param_ranges[j][0] + param_val_rescale[j]  * (param_ranges[j][1] - param_ranges[j][0]))
                self.params.append(Param(self.param_list, self.init_param_file, param_val))
                self.params[-1].param_error(self.mol_list, self.ref_data, self.log, self.processors, self.ref_rele)
                params_added += 1
                self.params_rescale.append(param_val_rescale)
            params_tried += 1
        print params_added, ' parameters added in ', params_tried, ' attempts'   

    # Generate random parameters near a specific point
    def gen_params_near(self, number, startpt, distmax, param_ranges):
        if len(self.params_rescale) < len(self.params):
            self.param_rescale()
        distcut = distmax * 1 
        print 'Generating ',number,' parameter sets within ', '%.4f'%distcut, ' of ', startpt, self.params[startpt].param_val

        params_added, params_tried = 0, 0
        while params_added < number and params_tried < number*100:
            param_val_rescale = []
            for j in range(0, len(self.param_list)):
                param_val_rescale.append((random.random()*1.0 - 0.5) * distcut)
            dist = 0.0
            for j in range(0, len(self.param_list)):
                dist += param_val_rescale[j]**2
            dist = math.sqrt(dist)
            if dist < distcut:
                param_val = []
                for j in range(0, len(self.param_list)):
                    param_val.append(self.params[startpt].param_val[j] + param_val_rescale[j]  * (param_ranges[j][1] - param_ranges[j][0]))
                self.params.append(Param(self.param_list, self.init_param_file, param_val))
                self.params[-1].param_error(self.mol_list, self.ref_data, self.log, self.processors, self.ref_rele)
                params_added += 1
                self.params_rescale.append(param_val_rescale)
            params_tried += 1
        print params_added, ' parameters added in ', params_tried, ' attempts'

    # Plot error
    def plot_error(self):
        error = []
        param_val = []
        for par in self.params:
            error.append(par.rms)
            param_val.append(par.param_val)
        error = np.array(error)
        param_val = np.array(param_val)

        fig = plt.figure()
        ax =fig.add_subplot(111,projection='3d')
        for i in range(0,len(error)):
            ax.scatter(param_val[i][0],param_val[i][1],error[i],marker='x',c='#ff0000')
        plt.show()

    # Plot error by step
    def plot_error_step(self):
        error = []
        i = []
        for par in self.params:
            i.append(len(i))
            error.append(par.rms)
        plt.scatter(i,error)
        plt.show()


    # Fit polynomial regression to parameter errors
    def fit_error(self):
        # Set up array of errors
        error = []
        param_val = []
        for par in self.params:
            error.append(par.rms)
            param_val.append(par.param_val)
        #print 'error', error
        #print 'param_val', param_val
        error = np.array(error)
        param_val = np.array(param_val)

        param_name = []
        for par in self.param_list:
            param_name.append(str(par[0])+'_'+par[1])


        plt.show()



class Param(object):
    def __init__(self, param_list, init_param_file, param_val):
        self.param_list = param_list
        self.init_param_file = init_param_file
        self.param_val = param_val

    def read_orbs(self, at_types):
        self.pfile = open(self.init_param_file, 'r')
        # Read the orbital parameters from INDO parameter file
        # orb_params contains atomic number, n, sp zeta, d zetas 1 & 2, and d coeffs 1 & 2
        self.params = []
        for i in range(0,len(at_types)):
            self.params.append([])

        pline = self.pfile.readline()
        for i in range(0,80):
            line = pline.split()
            if line[1] in at_types:
                # Read data for this atom
                k = at_types.index(line[1])
                self.params[k].append(int(line[0]))
                self.params[k].append(int(line[3]))
                pline = self.pfile.readline()
                line = pline.split()
                for j in range(0,5):
                    self.params[k].append(float(line[j]))
                for j in range(0,6):
                    pline = self.pfile.readline()

            else:
                for j in range(0,7):
                    pline = self.pfile.readline()


    # Edit external parameter file
    def edit_external(self):
        with open(self.init_param_file) as orig:
            olines = []
            for oline in orig:
                olines.append(oline)
        with open("external.txt",'w') as new:
            for i in range(len(olines)):
                new.write(olines[i])  
            for i in range(0,len(self.param_list)):
                new.write('     %-12sAg%17.6f\n'%(str(self.param_list[i][1]), self.param_val[i]))       




    # Run MOPAC
    def run_mopac(self, mol):
       # for mol in mol_list:
        os.system('/home/gieseking/Programs/MOPAC_clean/MOPAC2016.exe ' + mol)


    # Compute error for a reference set
    def param_error(self, mol_list, ref_data, log, processors, ref_rele):
       # self.edit_param()
        self.edit_external()
       # self.run_mopac(mol_list)
       # pool = mp.Pool(processes=4)
        pool = ThreadPool(processors)
        pool.map_async(self.run_mopac, [mol for mol in mol_list])
        pool.close()
        pool.join()

        refs = []

        log.write("Parameter values:\n")
        for par in self.param_list:
            par_string = str(par[0]) + ' ' + par[1]
            log.write(string.rjust(par_string,12))
        log.write('\n')
        s =''
        for val in self.param_val:
            log.write(string.rjust('%.6f'%val,12))
            s += string.rjust('%.6f'%val,10)
        log.write('\n')
        print s

        log.write("       Molecule    Ref E   Calc E    Error   Weight\n")
        #print ("       Molecule    Ref E   Calc E    Error   Weight")

        out_heat = []
        for i in range(0, len(mol_list)):
            prog = ' '
            # Read .out file & match states to reference data
            out, prog = open_mol(mol_list[i]+'.out', prog)
           # out.read_states(ref_data.ref_list[i].types)
           # ref_data.ref_list[i].match_refs(out.states)
            out.read_heat()  
            ref_data.ref_list[i].calc_error(out.heat)
            out_heat.append(out.heat)
        
           # # Append reference data to master ref_list
           # for exc in ref_data.ref_list[i].exc:
           #     r = [mol_list[i], exc.energy, exc.comp_energy, exc.error_energy, exc.weight, exc.match_found]
           #     refs.append(r)
           #     if r[5] == True:
           #         log.write("%15s %8.4f %8.4f %8.4f %8.4f\n"%(r[0],r[1],r[2],r[3],r[4]))
           #         #print ("%15s %8.4f %8.4f %8.4f %8.4f"%(r[0],r[1],r[2],r[3],r[4]))
           #     else:
           #         log.write("%15s %8.4f   **WARNING: No match found**\n"%(r[0],r[1]))
           #         #print ("%15s %8.4f   **WARNING: No match found**"%(r[0],r[1]))
            
           # Append reference data to master ref_list
            for hof in ref_data.ref_list[i].hof:
                r = [mol_list[i], hof.ref_heat, hof.comp_heat, hof.error_heat, hof.heat_weight, hof.match_found]
                refs.append(r)
                if r[5] == True:
                    log.write("%15s %8.4f %8.4f %8.4f %8.4f\n"%(r[0],r[1],r[2],r[3],r[4]))
                    #print ("%15s %8.4f %8.4f %8.4f %8.4f"%(r[0],r[1],r[2],r[3],r[4]))
                else:
                    log.write("%15s %8.4f   **WARNING: No match found**\n"%(r[0],r[1]))
    
        calc_rele = True 
        for eq_geo in ref_rele.eq_geos:
            if out_heat[eq_geo[2]] == 0:
                calc_rele = False
                break
            else:
                continue

        if calc_rele == True:          
            pm_rele = PM_Rel_E(out_heat, ref_rele)
            for i in range(0, len(ref_rele.ref_reles)):
               # if ref_rele.ref_reles[i][4] != 80.0:
                   # if (ref_rele.ref_reles[i][1] * pm_rele.pm_reles[i]) < 0:
                       # pm_rele.rele_error[i] = pm_rele.rele_error[i] * 100 
                r = [ref_rele.ref_reles[i][0], ref_rele.ref_reles[i][1], pm_rele.pm_reles[i], pm_rele.rele_error[i], ref_rele.ref_reles[i][4], calc_rele] 
                refs.append(r)
                if r[2] != 0:
                    log.write("%15s %8.4f %8.4f %8.4f %8.4f\n"%(r[0],r[1],r[2],r[3],r[4]))
                else:
                    log.write("%15s %8.4f   **WARNING: No match found**\n"%(r[0],r[1]))
 
        # Compute errors: Weighted MSE and Weighted RMS
        mse = 0.0
        mae = 0.0
        rms = 0.0
        count = 0
        wtcount = 0.0
        penalty = 40000.0
        
        for i in refs:
            if i[5] == True:
                mse += i[3]*i[4]
                mae += abs(i[3])*i[4]
                rms += i[3]**2*i[4]
                count += 1
                wtcount += i[4]
            else:
                mse += penalty*i[4]
                mae += penalty*i[4]
                rms += penalty**2*i[4]
                count += 1
                wtcount += i[4]
        
        mse = mse/wtcount
        mae = mae/wtcount
        rms = math.sqrt(rms/wtcount)
        
        log.write('\n')
        #print ''
        if count < len(refs):
            print ("WARNING: %4i references not found and excluded in statistics"%(len(refs)-count))
        #print ("Weighted mean signed error:      %8.4f" % mse)
        #print ("Weighted mean absolute error:    %8.4f" % mae)
        print ("              RMS error: %8.4f" % rms)
        log.write("Weighted mean signed error:      %8.4f\n" % mse)
        log.write("Weighted mean absolute error:    %8.4f\n" % mae)
        log.write("Weighted root mean square error: %8.4f\n\n" % rms)

        self.mse = mse
        self.mae = mae
        self.rms = rms

#        calc_rele = True
#        for eq_geo in ref_rele.eq_geos:
#            if out_heat[eq_geo[2]] == 0:
#                calc_rele = False
#                break
#            else:
#                continue
#        if calc_rele == True:
#           pm_rele = PM_Rel_E(out_heat, ref_rele)
#
#        log.write("       Molecule    Ref Rela E   Calc Rela E    Error   Weight\n")
#        rels = []
#        for i in range(0, len(ref_rele.ref_reles)):
#            r = [ref_rele.ref_reles[i][0], ref_rele.ref_reles[i][1], pm_rele.pm_reles[i], pm_rele.rele_error[i], pm_rele.weight]
#            rels.append(r)
#            log.write("%15s %8.4f %8.4f %8.4f %8.4f\n"%(r[0],r[1],r[2],r[3],r[4]))
#
#        mse = 0.0
#        mae = 0.0
#        rms = 0.0
#        count = 0
#        wtcount = 0.0
#
#        for i in rels:
#            mse += i[3]*i[4]
#            mae += abs(i[3])*i[4]
#            rms += i[3]**2*i[4]
#            count += 1
#            wtcount += i[4]
#
#        mse = mse/wtcount
#        mae = mae/wtcount
#        rms = math.sqrt(rms/wtcount)
#
#        log.write('\n')
#        #print ''
#        if count < len(refs):
#            print ("WARNING: %4i references not found and excluded in statistics"%(len(refs)-count))
#        #print ("Weighted mean signed error:      %8.4f" % mse)
#        #print ("Weighted mean absolute error:    %8.4f" % mae)
#        print ("              RMS error: %8.4f" % rms)
#        log.write("Weighted mean signed error:      %8.4f\n" % mse)
#        log.write("Weighted mean absolute error:    %8.4f\n" % mae)
#        log.write("Weighted root mean square error: %8.4f\n\n" % rms)






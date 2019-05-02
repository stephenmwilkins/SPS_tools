

import cPickle as pickle
import numpy as np
from scipy import ndimage

parameter_labels = {}
# parameter_labels['log10age'] = r'\log_{10}(t_{SF}/Myr)'
parameter_labels['log10n_H'] = r'\log_{10}(n_H/cm^{-2})'
parameter_labels['log10U'] = r'\log_{10}(U)'
parameter_labels['CO'] = r'\log_{10}[(C/O)/(C/O)_{\odot}]'
parameter_labels['d2m'] = r'\xi_{d}'




class grid:

    def __init__(self, SPS, IMF, lines, nebular_grid_dir = '../../BuildGrid/Full/grids/', stellar_grid_dir = '../../../stellar/BuildGrid/constant/grids/', generate_SED = False):

        # --- read in the nebular line grids
        
        model = '/'.join([SPS, IMF])

        self.lines = lines
                
        self.nebular_line_grid = False        
        
        
        print '--------------------- reading nebular line grid ---------------------'
        
        print '---- lines:'
        
        for line in self.lines:
        
            for l in line.split(','):
    
                print l
    
                iline_grid = pickle.load(open(nebular_grid_dir + model + '/LINES/'+l+'.p','r'))
    
                if not self.nebular_line_grid: self.nebular_line_grid = iline_grid
    
                self.nebular_line_grid[l] = {}
                self.nebular_line_grid[l]['luminosity'] = iline_grid['luminosity']
                self.nebular_line_grid[l]['lam'] = iline_grid['lam']

        print '---- parameter ranges:'
        
        for parameter in ['log10Z','log10n_H','log10U','d2m','CO']:
        
            print parameter, self.nebular_line_grid[parameter]

        if generate_SED:

            print '--------------------- reading nebular continuum grid ---------------------'

            self.nebular_continuum_grid = pickle.load(open(nebular_grid_dir + model + '/nebular_continuum.p','r'))
        
            self.lam = self.nebular_continuum_grid['lam']
        
            print '--------------------- reading stellar continuum grid ---------------------'

            self.stellar_grid = pickle.load(open(stellar_grid_dir + model + '/stellar.p','r'))
        
        
        
        
        

    def line_luminosity(self, p, line):
    
        # --- return the luminosity for a given set of parameters
        
        nebular_params = [[np.interp(p[parameter], self.nebular_line_grid[parameter], range(len(self.nebular_line_grid[parameter])))] for parameter in ['log10Z','log10n_H','log10U','d2m','CO']]
        
        intrinsic_luminosity = 1E7 * np.sum([10**ndimage.map_coordinates(self.nebular_line_grid[l]['luminosity'], nebular_params, order=1)[0] for l in line.split(',')])
        
        return intrinsic_luminosity
        
  
  
        
    def generate_nebular_continuum(self, p):  # --- generates nebular continuum SED
    
        iZ = np.interp(p['log10Z'], self.nebular_continuum_grid['log10Z'], range(len(self.nebular_continuum_grid['log10Z']))) # -- interpolated Z
            
        iZ1 = int(np.floor(iZ))
        iZ2 = int(np.ceil(iZ))
        
        w1 = 1. - (iZ - iZ1)
        w2 = 1. - (iZ2 - iZ) 
            
        if iZ1 == iZ2: w2 = 0.0 # necessary or they both get set to 1.0!
            
#         print iZ1, iZ2
#         print w1, w2
        
        self.intrinsic_nebular_continuum = 1E7 * (1.-p['fesc']) * (w1*self.nebular_continuum_grid['continuum'][iZ1, :]+w2*self.nebular_continuum_grid['continuum'][iZ2, :])
        
        # -- nebular emission is just calculated for 10 Myr constant SF hence 1E7 factor
  
 
  
    def generate_stellar(self, p):  # --- generates nebular continuum SED
    
        iZ = np.interp(p['log10Z'], self.nebular_continuum_grid['log10Z'], range(len(self.nebular_continuum_grid['log10Z']))) # -- interpolated Z
            
        iZ1 = int(np.floor(iZ))
        iZ2 = int(np.ceil(iZ))
        
        w1 = 1. - (iZ - iZ1)
        w2 = 1. - (iZ2 - iZ) 
        
        if iZ1 == iZ2: w2 = 0.0 # necessary or they both get set to 1.0!
        
        nage = (np.abs(self.stellar_grid['log10age'] - p['log10age'])).argmin() # -- nearest Z
        
        self.intrinsic_stellar = (w1*self.stellar_grid['L_nu'][nage, iZ1, :]+w2*self.stellar_grid['L_nu'][nage, iZ1, :])
      
  
    def EW(self, p, line):  # --- generates nebular continuum SED
    
    
        average_lam = np.mean([self.nebular_line_grid[l]['lam'] for l in line.split(',')])
    
        total_continuum = np.interp(average_lam, self.lam, self.intrinsic_stellar) + np.interp(average_lam, self.lam, self.intrinsic_nebular_continuum)  # erg s^-1 Hz^-1
      
        total_continuum *= 3E8 / (1E-10*average_lam**2) 
      
        EW = self.line_luminosity(p, line) / total_continuum
  
        return EW
  
  
  
  
        


import cPickle as pickle
import numpy as np
from scipy import ndimage
from scipy import stats

parameter_labels = {}
parameter_labels['fesc'] = r'f_{esc}'
parameter_labels['log10age'] = r'\log_{10}(t_{SF}/yr)'
parameter_labels['log10n_H'] = r'\log_{10}(n_H/cm^{-2})'
parameter_labels['log10U'] = r'\log_{10}(U_S)'
parameter_labels['CO'] = r'\log_{10}[(C/O)/(C/O)_{\odot}]'
parameter_labels['d2m'] = r'\xi_{d}'




class grid:

    def __init__(self, SPS, IMF, lines = False, nebular_grid_dir = '../../BuildGrid/grids/', stellar_grid_dir = '../../../stellar/BuildGrid/constant/grids/', generate_SED = False):

        # --- read in the nebular line grids
        
        model = '/'.join([SPS, IMF])

        self.lines = lines
                
        self.nebular_line_grid = False        
        
        
        print '--------------------- reading nebular line grid ---------------------'
        
        print '---- lines:'
        
        if lines:
        
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
        
        
        
        
        

    def generate_line_luminosity(self, p, line, dust_model = 'CF00'):
    
        # --- return the luminosity for a given set of parameters
        
        nebular_params = [[np.interp(p[parameter], self.nebular_line_grid[parameter], range(len(self.nebular_line_grid[parameter])))] for parameter in ['log10Z','log10n_H','log10U','d2m','CO']]
        
        self.intrinsic_line_luminosity = (1. - p['fesc']) * 1E7 * np.sum([10**ndimage.map_coordinates(self.nebular_line_grid[l]['luminosity'], nebular_params, order=1)[0] for l in line.split(',')])
        
        self.average_lam = np.mean([self.nebular_line_grid[l]['lam'] for l in line.split(',')])
        
        
        
        
        
        
        if dust_model == 'CF00':
        
            tau_BC = (1. - p['mu']) * p['tau_V_hat'] * (self.average_lam/5500.)**-1.3
    
            tau_ISM = p['mu'] * p['tau_V_hat'] * (self.average_lam/5500.)**-0.7
    
            T = np.exp(-( tau_BC + tau_ISM ))
        
            self.line_luminosity = self.intrinsic_line_luminosity * T


        elif dust_model == 'Calzetti':

            tau_nebular = p['tau_V_nebular'] * Calzetti(self.average_lam)/Calzetti(5500.) 
            T_nebular =  np.exp(-tau_nebular)

            self.line_luminosity = self.intrinsic_line_luminosity * T_nebular

        else:
        
            self.line_luminosity = self.intrinsic_line_luminosity
        
  
  
        
    def generate_nebular_continuum(self, p, dust_model = 'CF00'):  # --- generates nebular continuum SED
    
        iZ = np.interp(p['log10Z'], self.nebular_continuum_grid['log10Z'], range(len(self.nebular_continuum_grid['log10Z']))) # -- interpolated Z
            
        iZ1 = int(np.floor(iZ))
        iZ2 = int(np.ceil(iZ))
        
        w1 = 1. - (iZ - iZ1)
        w2 = 1. - (iZ2 - iZ) 
            
        if iZ1 == iZ2: w2 = 0.0 # necessary or they both get set to 1.0!
                    
        self.intrinsic_nebular_continuum = 1E7 * (1.-p['fesc']) * (w1*self.nebular_continuum_grid['continuum'][iZ1, :]+w2*self.nebular_continuum_grid['continuum'][iZ2, :])  # -- nebular emission is just calculated for 10 Myr constant SF hence 1E7 factor
         
        if dust_model == 'CF00':
        
            tau_BC = (1. - p['mu']) * p['tau_V_hat'] * (self.lam/5500.)**-1.3
    
            tau_ISM = p['mu'] * p['tau_V_hat'] * (self.lam/5500.)**-0.7
    
            T = np.exp(-( tau_BC + tau_ISM ))
        
        
            self.nebular_continuum = self.intrinsic_nebular_continuum * T
            
            
        elif dust_model == 'Calzetti':

            tau_nebular = p['tau_V_nebular'] * np.array([Calzetti(l)/Calzetti(5500.) for l in self.lam])
            T_nebular =  np.exp(-tau_nebular)

            self.nebular_continuum = self.intrinsic_nebular_continuum * T_nebular
        
        else:
        
            self.nebular_continuum = self.intrinsic_nebular_continuum
        
       
  
 
  
    def generate_stellar(self, p, dust_model = 'CF00'):  # --- generates nebular continuum SED
    
        iZ = np.interp(p['log10Z'], self.nebular_continuum_grid['log10Z'], range(len(self.nebular_continuum_grid['log10Z']))) # -- interpolated Z
            
        iZ1 = int(np.floor(iZ))
        iZ2 = int(np.ceil(iZ))
        
        w1 = 1. - (iZ - iZ1)
        w2 = 1. - (iZ2 - iZ) 
        
        if iZ1 == iZ2: w2 = 0.0 # necessary or they both get set to 1.0!
        
        nage = (np.abs(self.stellar_grid['log10age'] - p['log10age'])).argmin() # -- nearest Z
        
        self.intrinsic_stellar = (w1*self.stellar_grid['L_nu'][nage, iZ1, :]+w2*self.stellar_grid['L_nu'][nage, iZ2, :])
      
      
      
        if dust_model == 'CF00':
        
            intrinsic_stellar_BC = (w1*self.stellar_grid['L_nu'][0, iZ1, :]+w2*self.stellar_grid['L_nu'][0, iZ2, :])
        
            intrinsic_stellar_old = self.intrinsic_stellar - intrinsic_stellar_BC
            
            tau_BC = (1. - p['mu']) * p['tau_V_hat'] * (self.lam/5500.)**-1.3
    
            tau_ISM = p['mu'] * p['tau_V_hat'] * (self.lam/5500.)**-0.7
    
            T_BC = np.exp(-( tau_BC + tau_ISM ))
            T_ISM = np.exp(-( tau_ISM ))
        
            self.stellar = p['fesc']*self.intrinsic_stellar + (1.-p['fesc'])*(intrinsic_stellar_BC*T_BC + intrinsic_stellar_old*T_ISM)
        
        
        elif dust_model == 'Calzetti':

            tau_nebular = p['tau_V_nebular'] * np.array([Calzetti(l)/Calzetti(5500.) for l in self.lam])
            tau_stellar = 0.44 * tau_nebular
            T_stellar =  np.exp(-tau_stellar)

            self.stellar = p['fesc']*self.intrinsic_stellar + (1.-p['fesc'])*(self.intrinsic_stellar*T_stellar)
        
        else:
        
            self.stellar = self.intrinsic_stellar
      
      
    
    def generate_SED(self, p, dust_model = 'CF00'):
    
        
        self.generate_stellar(p, dust_model = dust_model)
        
        self.generate_nebular_continuum(p, dust_model = dust_model)
        
        # --- need to add lines
        
        self.intrinsic_total = self.intrinsic_stellar + self.intrinsic_nebular_continuum
    
        self.total = self.stellar + self.nebular_continuum
        
    
    
    
    def generate_UVCS(self, range = [1500.,3000.]):
    
        s = [(self.lam>range[0])&(self.lam<range[1])]
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(self.lam[s]),np.log10(self.intrinsic_total[s]))
        self.UVCS_int = slope - 2.
    
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(self.lam[s]),np.log10(self.total[s]))
        self.UVCS = slope - 2.
    
      
      
  
    def EW(self):  # --- generates nebular continuum SED
    
    
        total_continuum = np.interp(self.average_lam, self.lam, self.stellar) + np.interp(self.average_lam, self.lam, self.nebular_continuum)  # erg s^-1 Hz^-1
      
        total_continuum *= 3E8 / (1E-10*self.average_lam**2) 
      
        EW = self.line_luminosity / total_continuum
  
        return EW




    def fesc(self, l):
    
        # --- not correct - no line contribution - fine in the UV
    
        Fesc = self.total/self.intrinsic_total
        
        return np.interp(l, self.lam, Fesc)
        
  
    def A1500(self):
    
        return -2.5*np.log10(self.fesc(1500.))
                
 
 
 
  
def Calzetti(l):

    # --- Calzetti dust curve
    
    l_mu = l/1E4
    
    x = 1./l_mu 
     
    R_V_prime = 4.05
    
    
    if l_mu >= 0.63: k = 2.659 * (-1.857 + 1.040*x) + R_V_prime
    
    if l_mu < 0.63: k = 2.659 * (-2.156 + 1.509*x - 0.198*x**2 + 0.011*x**3) + R_V_prime
        
    return k
            
        
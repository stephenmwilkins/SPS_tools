

import numpy as np
import cPickle as pickle
from scipy import ndimage

age = 1E9 # --- duration of constant SF / yr 
n = 100000 # --- number of star particles
Z = 0.001 # --- metallicity

ages = age*np.random.random(n) # --- assume ages are uniformly spread over the age of the galaxy
metallicities = Z*np.ones(n) # --- assume they all have the same metallicity
mass = age/float(n) # --- mass of each star particle to give SFR = 1 Msol/yr 


line = 'HI6563'

grid = 'BPASSv2.1.binary_ModSalpeter_300'
# grid = 'BPASSv2.1.binary_ModSalpeter_100'
# grid = 'P2_ModSalpeter_100'

print '--------',line

line_grid = pickle.load(open('../'+grid+'/line_grid.p','r'))  # --- open grid

line_lam = line_grid[line]['lam']

line_luminosity = 0.0
stellar_continuum = 0.0
nebular_continuum = 0.0

for age, metallicity in zip(ages, metallicities):

    p = {'log10age': np.log10(age), 'log10Z': np.log10(metallicity)}

    params = [[np.interp(p[parameter], line_grid[parameter], range(len(line_grid[parameter])))] for parameter in ['log10age','log10Z']] # used in interpolation
    
    line_luminosity += mass * 10**ndimage.map_coordinates(line_grid[line]['luminosity'], params, order=1)[0]
    stellar_continuum += mass * ndimage.map_coordinates(line_grid[line]['stellar_continuum'], params, order=1)[0]
    nebular_continuum += mass * ndimage.map_coordinates(line_grid[line]['nebular_continuum'], params, order=1)[0]
    
print line_luminosity        
print stellar_continuum
print nebular_continuum

total_continuum = (stellar_continuum + nebular_continuum)*3E8/(line_lam**2 * 1E-10) 

EW = line_luminosity/total_continuum

print EW


import numpy as np
import cPickle as pickle
from scipy import ndimage

age = 1E8 # --- duration of constant SF / yr 
n = 100000 # --- number of star particles
Z = 0.001 # --- metallicity

ages = age*np.random.random(n) # --- assume ages are uniformly spread over the age of the galaxy
metallicities = Z*np.ones(n) # --- assume they all have the same metallicity
mass = age/float(n) # --- mass of each star particle to give SFR = 1 Msol/yr 

grid = 'BPASSv2.1.binary_ModSalpeter_300'
# grid = 'BPASSv2.1.binary_ModSalpeter_100'
# grid = 'P2_ModSalpeter_100'

f = 'FAKE.FAKE.1500'

print '--------',f

L_grid = pickle.load(open('../'+grid+'/FILTERS/'+f+'.p','r'))  # --- open grid

L = {}
for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: L[t] = 0.0

for age, metallicity in zip(ages, metallicities):

    p = {'log10age': np.log10(age), 'log10Z': np.log10(metallicity)}

    params = [[np.interp(p[parameter], L_grid[parameter], range(len(L_grid[parameter])))] for parameter in ['log10age','log10Z']] # used in interpolation
                 
    for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: L[t] += mass * ndimage.map_coordinates(L_grid[t], params, order=1)[0] # interpolate grid


for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: print t, L[t] # --- print luminosity from each component

print '-----'

print 'ratio of nebular to total:', L['nebular']/L['total'] # --- print ratio of nebular to total

# --- determine the approximate SFR calibration for the FUV. Based on Kennicutt & Evans 2012 this should be around 10^43.35 and you'd expect a value around this for 100 Myr constant SF.

print 'FUV SFR calibration:', np.log10(L['total']*3E8/(1500E-10)) 
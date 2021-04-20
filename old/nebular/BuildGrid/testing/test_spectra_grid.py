

import numpy as np
import cPickle as pickle
from scipy import ndimage
from scipy import integrate

age = 1E8
n = 1000
Z = 0.001


ages = age*np.random.random(n)
metallicities = Z*np.ones(n)
mass = age/float(n) # mass of each star particle to give SFR = 1 Msol/yr

grid = 'BPASSv2.1.binary_ModSalpeter_300'





SED_grid = pickle.load(open('../grids/'+grid+'/SED_grid.p','r'))

lam = SED_grid['lam']

SED_grid['log10age'] = 10**SED_grid['log10age']

L = {}
for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: L[t] = np.zeros(len(lam))


for age, metallicity in zip(ages, metallicities):

    p = {'log10age': np.log10(age), 'log10Z': np.log10(metallicity)}

    params = [[np.interp(p[parameter], SED_grid[parameter], range(len(SED_grid[parameter])))] for parameter in ['log10age','log10Z']]
             
    for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: L[t] += mass * np.array([ndimage.map_coordinates(SED_grid[t][:,:,i], params, order=1)[0] for i,l in enumerate(lam)])

    
    

filters = ['FAKE.FAKE.1500']

filters_dir = '/home/s/sw/sw376/Dropbox/Research/Utilities/filters/'


for f in filters:

    print '-----',f

    d = np.loadtxt(filters_dir+'/'.join(f.split('.'))+'.txt').T

    T = np.interp(lam, d[0], d[1]) # --- interpolate filter transmission grid on to redshifted SED grid 
    T[T<0.01] = 0.0

    for t in ['total', 'stellar', 'nebular', 'nebular_continuum']: 
        print t, integrate.trapz((1./lam) * L[t] * T, x = lam) / integrate.trapz((1./lam) * T, x = lam)

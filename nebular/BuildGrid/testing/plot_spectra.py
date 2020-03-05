

import numpy as np
import cPickle as pickle
from scipy import ndimage
from scipy import integrate

import matplotlib.pyplot as plt
import matplotlib.cm as cm



plt.style.use('simple')

fig = plt.figure(figsize = (4,4))

left  = 0.15 
bottom = 0.15 
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))



model = 'BPASSv2.1.binary/ModSalpeter_300'
Z = 0.001
log10age = 6.5

nebular_grid = pickle.load(open('../grids/'+model+'/nebular.p','r'))


ia = (np.abs(nebular_grid['log10age'] - log10age)).argmin()
iZ = (np.abs(nebular_grid['Z'] - Z)).argmin()


print ia, nebular_grid['log10age'][ia]
print iZ, nebular_grid['Z'][iZ]


ax.plot(np.log10(nebular_grid['lam']), np.log10(nebular_grid['stellar'][ia,iZ]), c='k', lw=4, alpha=0.2) 
ax.plot(np.log10(nebular_grid['lam']), np.log10(nebular_grid['stellar'][ia,iZ]+nebular_grid['nebular'][ia,iZ])) 

print np.max(nebular_grid['stellar'][ia,iZ])
print np.max(nebular_grid['nebular_continuum'][ia,iZ])



# --- compare with original stellar grid

stellar_grid = pickle.load(open('../../../../stellar/BuildGrid/SSP/grids/'+model+'/stellar.p','r'))

print stellar_grid.keys()

ia = (np.abs(stellar_grid['log10age'] - log10age)).argmin()
iZ = (np.abs(stellar_grid['Z'] - Z)).argmin()

print len(stellar_grid['lam']), np.min(stellar_grid['lam']), np.max(stellar_grid['lam'])


print ia, stellar_grid['log10age'][ia]
print iZ, stellar_grid['Z'][iZ]

ax.plot(np.log10(stellar_grid['lam']), np.log10(stellar_grid['L_nu'][ia,iZ]), c='k', lw=1) 



ax.set_xlabel(r'$\rm \log_{10}(\lambda/\mu m)$')
ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\ M_{\odot}^{-1})$')

fig.savefig('spectra.pdf')

fig.clf()
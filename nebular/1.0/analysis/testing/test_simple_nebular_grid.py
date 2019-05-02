

import cPickle as pickle
import numpy as np
from scipy import ndimage




line = 'CIII1909,CIII1907'

Z = 0.001
log10Z = np.log10(Z)

age = 8.


SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

grid_name = '_'.join([SPS, IMF])

nebular_grid = pickle.load(open('../BuildGrid/grids/' + grid_name + '/simple_nebular.p','r'))

line_lam = np.mean(np.array([nebular_grid[l]['lam'] for l in line.split(',')]))

print line_lam

luminosity = 1E7 * np.sum(np.array([10**np.interp(log10Z, nebular_grid['log10Z'], nebular_grid[l]['luminosity']) for l in line.split(',')]))
 
print luminosity

Lnu = np.array([np.interp(log10Z, nebular_grid['log10Z'], nebular_grid['SED.continuum'][:, i]) for i,l in enumerate(nebular_grid['SED.lam'])])

nebular_continuum = 1E7 * np.interp(line_lam, nebular_grid['SED.lam'], Lnu)

print nebular_continuum






# 
# # --- constant built up from bursts  --- THIS WORKS
# 
# lam = np.load('../../stellar/SSP/outputs/'+SPS+'/lam.npy')
# 
# tot_Lnu = np.zeros(lam.shape)
# totSF = 0.0    
# 
# # --- determine the first 10 Myr first
# 
# for log10age in np.arange(6.,age+0.1,0.1):
# 
#     Lnu = np.load('../../stellar/SSP/outputs/'+SPS+'/'+IMF+'/'+str(Z)+'_'+str(log10age)+'.npy') # W Hz^-1 M_sol^-1
#     Lnu *= 1E7 # erg s^-1 Hz^-1 
#     
#     SF = 10**log10age - totSF # assumes constant SF
#     totSF += SF
#     tot_Lnu += SF * Lnu
#     
# 
# stellar_continuum = np.interp(line_lam, lam, tot_Lnu)
# 
# print stellar_continuum




# # --- constant built up from bursts from grid  --- THIS WORKS

stellar_grid = pickle.load(open('../../stellar/BuildGrid/grids/' + grid_name + '/stellar.p','r'))

lam = stellar_grid['lam']

tot_Lnu = np.zeros(lam.shape)
totSF = 0.0    

for ib,log10age in enumerate(np.arange(7.,age+0.1,0.2)):

    # Lnu = stellar_grid['L_nu'][ib, 2, :]
    Lnu = np.array([np.interp(log10Z, stellar_grid['log10Z'], stellar_grid['L_nu'][ib, :, i]) for i,l in enumerate(lam)])
    
    SF = 10**log10age - totSF # assumes constant SF
    totSF += SF
    tot_Lnu += SF * Lnu
    
    
stellar_continuum = np.interp(line_lam, lam, tot_Lnu)

print stellar_continuum




# --- determine total continuum and EW


continuum = (stellar_continuum + nebular_continuum) * 3E8 / (line_lam**2 * 1E-10) # f_lam


EW = luminosity / continuum

print EW
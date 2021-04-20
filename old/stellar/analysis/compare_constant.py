

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
import cPickle as pickle



plt.style.use('simple')

fig = plt.figure( figsize=(6,4) )

left  = 0.2  
bottom = 0.2   
height = 0.7
width = 0.7

ax = fig.add_axes((left, bottom, width, height))



# SPS = 'P2/ModSalpeter_100'
SPS = 'BPASSv2.binary/ModSalpeter_100'
SPS = 'BPASSv2.1.binary/ModSalpeter_300'

Z = 0.004
age = 8.0


# 
# 
# # --- constant SFH as produced by original code
# 
# lam = np.load('../CONST/outputs/'+SPS.split('/')[0]+'/lam.npy')
# 
# Lnu = np.load('../CONST/outputs/'+SPS+'/'+str(Z)+'_'+str(age)+'.npy')
# Lnu *= 1E7 # erg s^-1 Hz^-1 
# 
# # if SPS == 'P2/ModSalpeter_100': Lnu /= 1E6
# 
# s = [(lam<10000)&(lam>100)]
# 
# ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = 'original code', c = 'k', alpha = 0.3)
# 


# --- constant built up from bursts

lam = np.load('../SSP/outputs/'+SPS.split('/')[0]+'/lam.npy')

s = [(lam<10000)&(lam>1000)]

tot_Lnu = np.zeros(lam.shape)
totSF = 0.0    

# --- determine the first 10 Myr first

for log10age in np.arange(6.,age+0.1,0.1):

    Lnu = np.load('../SSP/outputs/'+SPS+'/'+str(Z)+'_'+str(log10age)+'.npy') # W Hz^-1 M_sol^-1
    Lnu *= 1E7 # erg s^-1 Hz^-1 
    
    SF = 10**log10age - totSF # assumes constant SF
    totSF += SF
    tot_Lnu += SF * Lnu
    

ax.plot(np.log10(lam[s]), np.log10(tot_Lnu[s]), label = 'reconstructed', c = 'k', alpha = 0.1, lw =2)



# --- constant from grid



stellar_grid = pickle.load(open('../BuildGrid/constant/grids/' + SPS + '/stellar.p','r'))

print stellar_grid['L_nu'].shape

iZ = int(np.interp(np.log10(Z),stellar_grid['log10Z'], range(len(stellar_grid['log10Z'])))) # -- interpolated Z

iage = int(np.interp(age, stellar_grid['log10ages'], range(len(stellar_grid['log10ages'])))) # -- interpolated Z      

print iZ, iage      

lam = stellar_grid['lam']

s = [(lam<10000)&(lam>1000)]

Lnu = stellar_grid['L_nu'][iage, iZ, :]
    
print np.max(Lnu[s])

ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = 'stellar constant grid', c = 'k', alpha = 0.8, lw =0.5)


iage = int(np.interp(age, stellar_grid['log10ages'], range(len(stellar_grid['log10ages'])))) - 1

print iZ, iage      

lam = stellar_grid['lam']

s = [(lam<10000)&(lam>1000)]

Lnu = stellar_grid['L_nu'][iage, iZ, :]
    
print np.max(Lnu[s])

ax.plot(np.log10(lam[s]), np.log10(Lnu[s]), label = 'stellar constant grid', c = 'k', alpha = 0.8, lw =0.5)






ax.legend()

ax.set_ylabel(r'$L_{\nu}$')
ax.set_xlabel(r'$\log_{10}(\lambda/\AA)$')
    
fig.savefig('figures/constant.pdf')
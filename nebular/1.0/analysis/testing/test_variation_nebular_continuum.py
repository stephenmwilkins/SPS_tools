
import numpy as np


# ---- check how much variation there is in the nebular continuum

SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

Z = 0.001
log10age = 6.0
log10Z = np.log10(Z)
log10U = -2.0
log10n_H = 2.0
d2m = 0.3
CO = 0.0

model = '../RunGrid/Full/coutputs/C17.00/'+SPS+'/'+IMF+'/'+'_'.join(map(str,[Z,log10age,log10U,log10n_H,d2m,CO]))
lam, nebular, linecont = np.loadtxt(model+'.cont', delimiter='\t', usecols = (0,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
nebular_continuum = nebular - linecont


default = {}

for l in [1500.,2500,3500.]: default[l] = np.interp(l, lam[::-1], nebular_continuum[::-1])


# --- test effect of metallicity

print '------- Z'

for Z in [0.00001,0.0001, 0.001, 0.004, 0.01]: 

    model = '../RunGrid/Full/coutputs/C17.00/'+SPS+'/'+IMF+'/'+'_'.join(map(str,[Z,log10age,log10U,log10n_H,d2m,CO]))
    lam, nebular, linecont = np.loadtxt(model+'.cont', delimiter='\t', usecols = (0,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
    nebular_continuum = nebular - linecont
   
    print Z,
    for l in [1500.,2500,3500.]: print np.interp(l, lam[::-1], nebular_continuum[::-1])/default[l],
    print
    


    
# --- test effect of log10U

print '------- log10U'

Z = 0.001
log10age = 6.0
log10Z = np.log10(Z)
log10U = -2.0
log10n_H = 2.0
d2m = 0.3
CO = 0.0

for log10U in [-4.,-3.,-2.,-1.]: 

    model = '../RunGrid/Full/coutputs/C17.00/'+SPS+'/'+IMF+'/'+'_'.join(map(str,[Z,log10age,log10U,log10n_H,d2m,CO]))
    lam, nebular, linecont = np.loadtxt(model+'.cont', delimiter='\t', usecols = (0,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
    nebular_continuum = nebular - linecont
   
    print log10U,
    for l in [1500.,2500,3500.]: print np.interp(l, lam[::-1], nebular_continuum[::-1])/default[l],
    print
    

    
# --- test effect of log10n_H

print '------- log10n_H'

Z = 0.001
log10age = 6.0
log10Z = np.log10(Z)
log10U = -2.0
log10n_H = 2.0
d2m = 0.3
CO = 0.0

for log10n_H in [2.0,3.0,4.0]: 

    model = '../RunGrid/Full/coutputs/C17.00/'+SPS+'/'+IMF+'/'+'_'.join(map(str,[Z,log10age,log10U,log10n_H,d2m,CO]))
    lam, nebular, linecont = np.loadtxt(model+'.cont', delimiter='\t', usecols = (0,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
    nebular_continuum = nebular - linecont
   
    print log10n_H,
    for l in [1500.,2500,3500.]: print np.interp(l, lam[::-1], nebular_continuum[::-1])/default[l],
    print
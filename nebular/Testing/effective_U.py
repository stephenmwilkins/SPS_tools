

import numpy as np
import os
from scipy import integrate
import math



Z = 0.01
SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

log10U_S_0 = -2.0
log10n_H = 2.5


input_dir = '../../stellar/SSP/outputs'
lam = np.load(input_dir+'/'+SPS+'/lam.npy') # wavelength in \AA
nu = 3E8/(lam*1E-10)  # frequency in Hz
nu_log10 = np.log10(nu)
h = 6.626070040E-34 # J s
c = 3E8 # m s-1
nu_limit = c/(91.2E-9)


def get_log10Q(Z, log10age):

    model_file = '_'.join([str(Z), str(log10age)])
    
    Lnu = np.load(input_dir+'/'+SPS+'/'+IMF+'/'+model_file+'.npy') # (erg s^-1 Hz^-1) / M_sol yr^-1
    Lnu_log10 = np.log10(Lnu+1E-99)
    Lnu_log10 -= np.max(Lnu_log10) # normalises to the maximum. This is OK because the actually luminosity input is just log10Q.
    Lnu_log10[(lam<10.)] = -20.

    f = lambda x: np.interp(x, nu[::-1], Lnu[::-1])/(h*x)
    return np.log10(integrate.quad(f, nu_limit, nu_limit*1000)[0]) # 1 M_sol mass star cluster

    

log10Q_6 = get_log10Q(Z, 6.0)

totSF = 0.0
totQ = 0.0
culmulative_U_S = 0.0

for log10age in np.arange(6.0,9.1,0.1):

    log10Q = get_log10Q(Z, log10age)

    log10U_S = log10U_S_0 + 0.33*(log10Q-log10Q_6)
    
    SF = 10**log10age - totSF # assumes constant SF
    totSF += SF
    
    totQ += SF*10**log10Q
    
    culmulative_U_S += SF*10**log10Q * 10**log10U_S
    
    
    print log10age, log10Q, log10U_S, SF, totQ, np.log10(culmulative_U_S/totQ)
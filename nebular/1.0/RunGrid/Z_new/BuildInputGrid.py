



import numpy as np
import utilities as u
import os
from astropy.io import ascii

cloudy_version = 'C17.00'


SPS = 'BPASSv2.2.1.binary'
IMFs = ['1p0_100','1p0_300','ModSalpeter_100','1p7_100','1p7_300','Salpeter_100','Chabrier_100','Chabrier_300']
IMFs = ['ModSalpeter_300']



if SPS == 'P2': Zs = [0.0004, 0.004, 0.008, 0.02, 0.05]  
if SPS in ['BPASSv2.1.binary','BPASSv2.1.single','BPASSv2.2.binary','BPASSv2.2.single','BPASSv2.2.1.binary']: Zs = [0.00001, 0.0001, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04]

log10ages = np.arange(6.,10.1,0.1)

log10n_H = 2.5
log10U_S_0 = -2.0
CO = 0.0
d2m = 0.3

dry_run = False

if not dry_run:
    if not os.path.exists('coutputs/'+cloudy_version+'/'+SPS): os.mkdir('coutputs/'+cloudy_version+'/'+SPS)

i = 1 # model ID

for IMF in IMFs:

    if not dry_run:

        if not os.path.exists('coutputs/'+cloudy_version+'/'+SPS+'/'+IMF): os.mkdir('coutputs/'+cloudy_version+'/'+SPS+'/'+IMF)

    # --- read in wavelength and determine frequency


    # --- define a "zero" age, reference metallicity log10Q
    
    CLOUDY_SED, log10Q_orig = u.create_CLOUDY_SED(SPS, IMF, 0.02, 6.)

    log10Q_6_solar = log10Q_orig


    summary = {'Z': [], 'log10age':[], 'log10U_S': [], 'log10Q': [], 'log10Q_orig': []}


    for Z in Zs:

        for log10age in log10ages:

            CLOUDY_SED, log10Q_orig = u.create_CLOUDY_SED(SPS, IMF, Z, log10age)
        
            log10U_S = log10U_S_0 + (log10Q_orig - log10Q_6_solar)/3.

            log10Q = u.determine_log10Q(log10U_S, log10n_H) # --- this is the input         
        
            print(Z, log10age, '|', log10U_S, log10Q, log10Q_orig )  
          
            summary['Z'].append(Z)
            summary['log10age'].append(log10age)
            summary['log10U_S'].append(log10U_S)
            summary['log10Q'].append(log10Q)
            summary['log10Q_orig'].append(log10Q_orig)
          
            if not dry_run: u.make_CLOUDY_input(CLOUDY_SED, log10Q, SPS, IMF, Z, log10age, log10U_S_0, log10n_H, CO, d2m, i = i, run = False, cloudy_version=cloudy_version)

            i += 1

    ascii.write(summary,'input_summary_'+IMF+'.dat', names=['Z', 'log10age', 'log10U_S', 'log10Q', 'log10Q_orig'])                                                              
    print('created '+str(i-1)+' models')


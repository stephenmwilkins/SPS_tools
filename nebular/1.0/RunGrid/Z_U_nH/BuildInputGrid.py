



import numpy as np
import utilities as u
import os

cloudy_version = 'C17.00'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_300'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_100'

# SPS = 'P2'
# IMF = 'ModSalpeter_100'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_300'

SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'



if not os.path.exists('coutputs/'+cloudy_version+'/'+SPS): os.mkdir('coutputs/'+cloudy_version+'/'+SPS)
if not os.path.exists('coutputs/'+cloudy_version+'/'+SPS+'/'+IMF): os.mkdir('coutputs/'+cloudy_version+'/'+SPS+'/'+IMF)





    
if SPS == 'P2': Zs = [0.0004, 0.004, 0.008, 0.02, 0.05]  
if SPS in ['BPASSv2.1.binary','BPASSv2.1.single','BPASSv2.2.binary','BPASSv2.2.single']: Zs = [0.00001, 0.0001, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04]


# Zs = [0.004]

log10ages = np.array([6.0])

log10n_Hs = np.array([1.0,1.5,2.0,2.5,3.0,3.5,4.0])
log10U_S_0s =  np.array([-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0])
CO = 0.0
d2m = 0.3


dry_run = False

i = 1 # model ID

# --- read in wavelength and determine frequency

for Z in Zs:

    for log10age in log10ages:

        CLOUDY_SED, log10Q_orig = u.create_CLOUDY_SED(SPS, IMF, Z, log10age)

        if log10age == 6.0: log10Q_6 = log10Q_orig
        
        for log10U_S_0 in log10U_S_0s:
        
            log10U_S = log10U_S_0 + (log10Q_orig - log10Q_6)/3.

            for log10n_H in log10n_Hs:

                log10Q = u.determine_log10Q(log10U_S, log10n_H) # --- this is the input         
        
                print Z, log10age, '|', log10U_S, log10Q   
          
                if not dry_run: u.make_CLOUDY_input(CLOUDY_SED, log10Q, SPS, IMF, Z, log10age, log10U_S_0, log10n_H, CO, d2m, i = i, run = False, cloudy_version=cloudy_version)

                i += 1


                            
                            
                                
print 'created '+str(i-1)+' models'






import numpy as np
import utilities as u
import os
from shutil import copyfile

cloudy_version = 'C17.00'

SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_100'

# SPS = 'P2'
# IMF = 'ModSalpeter_100'


    
if SPS == 'P2': Zs = [0.0004, 0.004, 0.008, 0.02, 0.05]  
if SPS in ['BPASSv2.1.binary','BPASSv2.1.single']: Zs = [0.00001, 0.0001, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04]

log10ages = np.arange(6.,8.1,0.1)

log10n_Hs = np.arange(1.5,3.51,1.0)
log10U_S_0s = np.arange(-4.0,-0.5,1.0)





i = 1 # model ID

zero = 0

# --- read in wavelength and determine frequency

for Z in Zs:

    for log10age in log10ages:
        
        for log10U_S_0 in log10U_S_0s:        
            
            for log10n_H in log10n_Hs:

                try:

                    outdir = 'coutputs/'
                    output_file = outdir+cloudy_version+'/'+SPS+'/'+IMF+'/'+'_'.join([str(Z), str(log10age), str(log10U_S_0), str(log10n_H)])

                    statinfo = os.stat(output_file+'.lines')

                    if statinfo.st_size<1000: 
                        print i, statinfo.st_size, '|', Z, log10age, log10n_H, log10U_S_0
                        zero += 1
                        # os.system('qsub -t '+str(i)+' RunGrid.job')


                except:
                
                    print 'not found', i

                prevfile = output_file

                i += 1


                            
                            
                                
print 'created '+str(i-1)+' models'
print 'models that failed: '+str(zero)


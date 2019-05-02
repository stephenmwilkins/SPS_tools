



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

try:
    os.mkdir('coutputs/'+cloudy_version+'/'+SPS)
except:
    print 'already made SPS dir'

try:
    os.mkdir('coutputs/'+cloudy_version+'/'+SPS+'/'+IMF)
except:
    print 'already made IMF dir'


    
if SPS == 'P2': Zs = [0.0004, 0.004, 0.008, 0.02, 0.05]  
if SPS in ['BPASSv2.1.binary','BPASSv2.1.single']: Zs = [0.00001, 0.0001, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04]

log10ages = np.arange(6.,7.1,0.1)

log10n_Hs = np.arange(2.0,4.1,1.0)
log10Us = np.arange(-4.0,-0.9,0.5) 
d2ms = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
COs = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25]




i = 1 # model ID

zero = 0

# --- read in wavelength and determine frequency

for Z in Zs:

    for log10age in log10ages:
                
        for log10n_H in log10n_Hs:

            for log10U in log10Us:
      
                for d2m in d2ms:
        
                    for CO in COs:

                        try:

                            outdir = '/lustre/scratch/astro/sw376/SED/nebular/RunGrid/Full_NewZ/coutputs/'
                            output_file = outdir+cloudy_version+'/'+SPS+'/'+IMF+'/'+'_'.join([str(Z), str(log10age), str(log10U), str(log10n_H), str(d2m), str(CO)])

                            statinfo = os.stat(output_file+'.lines')

                            if statinfo.st_size<100: 
                                print i, statinfo.st_size, '|', Z, log10age, log10n_H, log10U, d2m, CO
                                zero += 1
            
#                                 copyfile(prevfile+'.cont', output_file+'.cont')
#                                 copyfile(prevfile+'.lines', output_file+'.lines')

    
                        except:
                        
                            print 'not found', i

                        prevfile = output_file

                        i += 1


                            
                            
                                
print 'created '+str(i-1)+' models'
print 'models that failed: '+str(zero)


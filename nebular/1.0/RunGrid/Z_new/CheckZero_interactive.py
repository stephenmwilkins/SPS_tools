



import numpy as np
import utilities as u
import os
from shutil import copyfile

cloudy_version = 'C17.00'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_300'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_100'

# SPS = 'P2'
# IMF = 'ModSalpeter_100'

SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'

    
if SPS == 'P2': Zs = [0.0004, 0.004, 0.008, 0.02, 0.05]  
if SPS in ['BPASSv2.1.binary','BPASSv2.1.single','BPASSv2.2.binary','BPASSv2.2.single']: Zs = [0.00001, 0.0001, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.03, 0.04]

log10ages = np.arange(6.,9.1,0.1)



outdir = 'coutputs/'+cloudy_version+'/'+SPS+'/'+IMF+'/'


i = 1 # model ID

failed = 0

# --- read in wavelength and determine frequency


for Z in Zs:

    for log10age in log10ages:
          
        output_file = outdir+'_'.join([str(Z), str(log10age)])

        if os.path.exists(output_file+'.lines'):

            statinfo = os.stat(output_file+'.lines')

            if statinfo.st_size<1000: 
                print i, statinfo.st_size, '|', Z, log10age
                failed += 1
                # os.system('/home/s/sw/sw376/software/Cloudy/c17.00/source/cloudy.exe < cinputs/'+str(i)+'.in')

                os.system('cp '+prev_file+'.lines '+output_file+'.lines')
                os.system('cp '+prev_file+'.cont '+output_file+'.cont')


        else:
        
            print 'not found', i, Z, log10age
            # os.system('/home/s/sw/sw376/software/Cloudy/c17.00/source/cloudy.exe < cinputs/'+str(i)+'.in')
            
            os.system('cp '+prev_file+'.lines '+output_file+'.lines')
            os.system('cp '+prev_file+'.cont '+output_file+'.cont')
            
            failed += 1

        i += 1

        prev_file = output_file
                            
                            
                                
print 'created '+str(i-1)+' models'
print 'models that failed: '+str(failed)


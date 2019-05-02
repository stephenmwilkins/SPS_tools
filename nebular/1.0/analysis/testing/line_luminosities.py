import numpy as np
# import matplotlib.pyplot as plt
# from scipy import integrate
import readCLOUDY



SPS = 'BPASSv2.binary'
IMF = 'Salpeter'
Z = 0.004


line = 'HI6563'
line = 'CIII1910'
# line = 'CIII1907'

totL = 0.0
tot_continuum = 0.0

totSF = 0.0

for log10age in np.arange(6,8.1,0.2):

    SF = 10**log10age - totSF
    
    totSF += SF
    
    model = '_'.join([IMF, str(Z), str(log10age)])
    
    lam, ID, emergent, continuum_nebular, continuum_stellar = readCLOUDY.nebular(SPS, model, lines = [line])
    
    totL += SF * (10**emergent[(ID==line)][0])

    tot_continuum += SF * (10**continuum_stellar[(ID==line)][0]+10**continuum_nebular[(ID==line)][0])

    EW = lam[(ID==line)][0]*totL/tot_continuum

    continuum_stellar2 = readCLOUDY.stellar(SPS, model, ID, lam)
    
    print log10age, SF, totSF, np.log10(totL)-6., continuum_stellar[0], np.log10(continuum_stellar2[0]), EW

    



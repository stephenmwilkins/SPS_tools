
import os
import sys
import numpy as np
import pickle
import read
from scipy import stats
from scipy import integrate



h = 6.626E-34	
c = 3.E8

cloudy_grid = 'Z_U_nH'

SPS = 'BPASSv2.2.binary'
IMF = 'ModSalpeter_300'

# SPS = 'BPASSv2.1.binary'
# IMF = 'ModSalpeter_100'
# 
# SPS = 'P2'
# IMF = 'ModSalpeter_100'
   
if not os.path.isdir('grids/'+SPS): os.mkdir('grids/'+SPS)
if not os.path.isdir('grids/'+SPS+'/'+IMF): os.mkdir('grids/'+SPS+'/'+IMF)    


# --- open input stellar grid (using for wavelength grid and rescaling)

stellar_grid = pickle.load(open('../../../stellar/BuildGrid/grids/'+SPS+'/'+IMF+'/stellar.p','rb'))

lam = stellar_grid['lam']



# --- determine log10Q_stellar to normalise

# based on default model with Z=0.004 and all lines with R>-2.0

lines = []
lines += ['MgII2803','MgII2796']
lines += ['ArIII7136','ArIII7751','ArIV4711','ArIV4740']
lines += ['NeIII3869','NeIII3967']
lines += ['FeIV3095','FeIV2836','FeIV2829']
lines += ['CII2325','CII2327','CIII1909','CIII1907','CIV1551','CIV1548','C-1335']
lines += ['OII2470','OII3729','OII3726','OIII4959','OIII5007','OIII2321','OIII4363','OIII1661','OIII1666']
lines += ['SiIII1892','SiIII1883','SiIII1206','SiIV1394','SiIV1403']
lines += ['NII6583']
lines += ['SII6731','SII6716','SIII9069','SIII9531','SIII3722','SIII6312']
lines += ['HeII1640','HeI10830','HeI3889','HeI3188','HeI2945','HeI2829','HeI20581','HeI5016','HeI3965','HeI7065','HeI5876','HeI4471','HeI4026','HeI3820','HeI3705','HeI6678','HeI4922','HeI18685']
lines += ['HI1216','HI1026','HI6563','HI4861','HI4340','HI4102','HI3970','HI3889','HI3835','HI3798','HI3771','HI3750','HI3734','HI3722','HI3712','HI3704','HI3697','HI3692','HI3687','HI3683','HI3679','HI3671','HI3669','HI18751','HI12818','HI10938','HI10049','HI9546','HI9229','HI9015','HI8863','HI8750','HI8665','HI8323','HI26251','HI21655','HI19445','HI18174']


lines = np.array(lines)


Zs = stellar_grid['Z']


log10n_Hs = np.arange(2.0,3.51,1.0)
log10U_S_0s = np.arange(-4.0,-0.5,1.0)


# log10n_Hs = np.array([2.0])
# log10U_S_0s = np.array([-2.0])

shape = (len(Zs), len(log10U_S_0s), len(log10n_Hs))

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# ---- get line fluxes and nebular contribution to continuum


for tmax in [7.0, 8.0]:

    line_grid = {}
    line_grid['Z'] = Zs
    line_grid['log10Z'] = np.log10(line_grid['Z'])
    line_grid['log10n_H'] = log10n_Hs
    line_grid['log10U_S_0'] = log10U_S_0s
    line_grid['lines'] = lines # --- list of lines

    for line in lines: 
    
        line_grid[line] = {}
        
        line_grid[line]['luminosity'] = np.zeros(shape)
        line_grid[line]['nebular_continuum'] = np.zeros(shape)
        line_grid[line]['stellar_continuum'] = np.zeros(shape)
        line_grid[line]['total_continuum'] = np.zeros(shape)

    for iZ, Z in enumerate(Zs):
    
        for iU, log10U_S_0 in enumerate(log10U_S_0s):
        
            for inH, log10n_H, in enumerate(log10n_Hs):
    

                print('-----------', Z, log10U_S_0, log10n_H)
    
                totSF = 0.0
    
                for log10age in np.arange(6.0, tmax+0.1, 0.1):

                    SF = 10**(np.min((log10age+0.05, tmax))) - totSF # --- SF associated with this age bin assumes constant SF
                    totSF += SF

                    

                    model = '_'.join([str(Z), str(log10age), str(log10U_S_0), str(log10n_H)])  # --- edit this for more complicated parameters

                    # ----- continuum

                    continuum_lam, continuum_nu, stellar, nebular_continuum, total, linecont = read.continuum(SPS, IMF, model, cloudy_grid = cloudy_grid, cloudy_version = 'C17.00')  # --- all L_nu

                    # ----- lines

                    line_lam, line_luminosity  = read.lines(SPS, IMF, model, lines = lines, cloudy_grid = cloudy_grid, cloudy_version = 'C17.00')

                    for i,line in enumerate(lines): 
    
                        line_grid[line]['lam'] = line_lam[i]
                        line_grid[line]['luminosity'][iZ,iU,inH] += SF*10**(line_luminosity[i]) 
                        line_grid[line]['nebular_continuum'][iZ,iU,inH] += SF*np.interp(line_lam[i], continuum_lam[::-1], nebular_continuum[::-1])
                        line_grid[line]['stellar_continuum'][iZ,iU,inH] += SF*np.interp(line_lam[i], continuum_lam[::-1], stellar[::-1])
                        line_grid[line]['total_continuum'][iZ,iU,inH] += SF*line_grid[line]['nebular_continuum'][iZ,iU,inH] + line_grid[line]['stellar_continuum'][iZ,iU,inH]

    
    
    tmax_label = {7.0:'10',8.0:'100'}

    pickle.dump(line_grid, open('grids/'+SPS+'/'+IMF+'/lines'+tmax_label[tmax]+'.p', 'w')) # --- save combined grid (can delete if necessary)


                            


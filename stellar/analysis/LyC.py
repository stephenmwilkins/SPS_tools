import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
import pickle


h = 6.626E-34	
c = 3.E8

plt.style.use('simple')

fig = plt.figure(figsize=(3.5,3.5))

left  = 0.15  
bottom = 0.15   
height = 0.8
width = 0.8


ax = fig.add_axes((left, bottom, width, height))






SPS = 'BPASSv2.2.1.binary/ModSalpeter_300'

grid = pickle.load(open('../BuildGrid/grids/'+SPS+'/stellar.p','rb'))


SPS_label = r'$\rm '+SPS+'$'


lam = grid['lam']

print(lam)

for iZ, Z, ls in zip(range(len(grid['Z'])), grid['Z'], ['-','-.','--',':']*3):

    color = cm.viridis(float(iZ)/float(len(grid['Z'])+1))

    N_LyC = []

    for ia,log10age in enumerate(grid['log10age']):

        # Lnu = np.load('../SSP/outputs/'+SPS+'/'+str(Z)+'_'+str(log10age)+'.npy') # erg s^-1 Hz^-1 M_sol
        
        Lnu = 1E-7 * grid['L_nu'][ia, iZ] # W s^-1 Hz^-1

        Llam = Lnu * c / (lam**2*1E-10) # W s^-1 \AA^-1
        
        nlam = (Llam*lam*1E-10)/(h*c) # s^-1 \AA^-1

        # nlam = (Llam*912*1E-10)/(h*c) # s^-1 \AA^-1

        f = lambda l: np.interp(l, lam, nlam)
        
        n_LyC = integrate.quad(f, 10.0, 912.0)[0]
            
#         nlam[lam>912.] = 0.0
# 
#         n_LyC_alt = integrate.trapz(nlam, x = lam)
        
        print(Z, log10age, np.log10(n_LyC))

        N_LyC.append(np.log10(n_LyC))

    
    if iZ==0: 
        Z_label = '0.00001' 
    else:
        Z_label = str(Z)
    
    ax.plot(grid['log10age'], np.array(N_LyC), c=color, lw=1, ls=ls, label=r'$\rm Z='+Z_label+'$')
    
    

ax.set_xlim([6.,9.])
ax.set_ylim([40,48.])

ax.text(0.5, 1.02, SPS_label, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, fontsize = 6, alpha = 0.5)




# ax.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5), fontsize = 6)

ax.legend(loc = 'upper right', fontsize = 6)




# left  = 0.25 
# bottom = 0.25   
# height = 0.3
# width = 0.3
# 
# total_ax = fig.add_axes((left, bottom, width, height))
# 
# total_ax.xaxis.set_tick_params(labelsize=4)
# total_ax.yaxis.set_tick_params(labelsize=4)
# 
# # ax.set_xlim([6.,9.])
# # ax.set_ylim([40,48.])















ax.set_xlabel(r'$\rm\log_{10}(age/yr)$')
ax.set_ylabel(r'$\rm\log_{10}(\dot{n}_{LyC}/s^{-1}M_{\odot}^{-1})$')
 
    
fig.savefig('figures/LyC.pdf')
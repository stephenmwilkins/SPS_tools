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


SPS_label = r'$\rm BPASSv2.2.1\ Modified\ Salpeter\ IMF\ m_{up}=300\ M_{\odot}$'


lam = grid['lam']


for t_max, ls, ms in zip([7.0,8.0],['-','--'],['o','^']):


    tot = []

    for iZ, Z in zip(range(len(grid['Z'])), grid['Z']):

        color = cm.viridis(float(iZ)/float(len(grid['Z'])+1))

        prevt = 0.0

        N_LyC = 0.0

        for ia,log10age in enumerate(np.arange(6.0, t_max+0.1, 0.1)):


            dt = 10**(np.min((log10age+0.05, t_max))) - prevt 
        
            prevt += dt

            # Lnu = np.load('../SSP/outputs/'+SPS+'/'+str(Z)+'_'+str(log10age)+'.npy') # erg s^-1 Hz^-1 M_sol
        
            Lnu = 1E-7 * grid['L_nu'][ia, iZ] # W s^-1 Hz^-1

            Llam = Lnu * c / (lam**2*1E-10) # W s^-1 \AA^-1
        
            nlam = (Llam*lam*1E-10)/(h*c) # s^-1 \AA^-1

            # nlam = (Llam*912*1E-10)/(h*c) # s^-1 \AA^-1

            f = lambda l: np.interp(l, lam, nlam)
        
            n_LyC = integrate.quad(f, 10.0, 912.0)[0]
             
            N_LyC += n_LyC * dt * (365.*24*60*60)
            
            # print Z, log10age, dt, prevt, np.log10(N_LyC)#
        
        tot.append(N_LyC)
    
        ax.scatter(np.log10(Z), np.log10(N_LyC), c=color, zorder =1, marker=ms, s=15)
    
    ax.plot(grid['log10Z'], np.log10(np.array(tot)), c='k', alpha=0.2, ls=ls, zorder =0, label=r'$\rm\log_{10}(t_{max})='+str(t_max)+'$')
    

ax.set_xlim([-5.1,-1.3])
# ax.set_ylim([40,48.])

ax.text(0.5, 1.02, SPS_label, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, fontsize = 6, alpha = 0.5)


ax.legend(fontsize=7)












ax.set_xlabel(r'$\rm\log_{10}(Z)$')
ax.set_ylabel(r'$\rm\log_{10}\left[\int_{0}^{t_{max}}\,\dot{n}_{LyC}(t)\,dt/M_{\odot}^{-1}\right]$')
 
    
fig.savefig('figures/LyC_total.pdf')
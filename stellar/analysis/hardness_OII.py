import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
import pickle


h = 6.626E-34	
c = 3.E8

h_eV = 4.135667662E-15

def lam_energy(E_eV):
    return (h_eV*c/E_eV)*1E10
    

print(lam_energy(13.6))
print(lam_energy(35.11))


plt.style.use('simple')

fig = plt.figure(figsize=(3.5,3.5))

left  = 0.15  
bottom = 0.15   
height = 0.8
width = 0.8


ax = fig.add_axes((left, bottom, width, height))


all_Z = np.array([0.00001, 0.0001, 0.001,0.002,0.004,0.006,0.008,0.01,0.02,0.03,0.04])


all_colors = np.array([cm.viridis(float(iZ)/float(len(all_Z)+1)) for iZ in range(len(all_Z))])

Zs = [0.004, 0.02]

SPS_models = ['BPASSv2.2.1.binary/ModSalpeter_300','BPASSv2.2.1.binary/ModSalpeter_100','P2/ModSalpeter_100']
line_styles = ['-','-.',':']

SPS_labels = {}
SPS_labels['BPASSv2.2.1.binary/ModSalpeter_300'] = r'$\rm BPASSv2.2.1\ Modified\ Salpeter\ IMF\ m_{up}=300\ M_{\odot}$'
SPS_labels['BPASSv2.2.1.binary/ModSalpeter_100'] = r'$\rm BPASSv2.2.1\ Modified\ Salpeter\ IMF\ m_{up}=100\ M_{\odot}$'
SPS_labels['P2/ModSalpeter_100'] = r'$\rm PEGASE.2\ Modified\ Salpeter\ IMF\ m_{up}=100\ M_{\odot}$'



for SPS, ls in zip(SPS_models,line_styles):

    print(SPS)

    grid = pickle.load(open('../BuildGrid/grids/'+SPS+'/stellar.p','rb'))
    
    lam = grid['lam']


    for Z in Zs:

        iZ = (np.abs(grid['Z'] - Z)).argmin()

        print(iZ, Z)

        color = all_colors[all_Z==Z][0]

        R = []

        for ia,log10age in enumerate(grid['log10age']):

            # Lnu = np.load('../SSP/outputs/'+SPS+'/'+str(Z)+'_'+str(log10age)+'.npy') # erg s^-1 Hz^-1 M_sol
        
            Lnu = 1E-7 * grid['L_nu'][ia, iZ] # W s^-1 Hz^-1

            Llam = Lnu * c / (lam**2*1E-10) # W s^-1 \AA^-1
        
            nlam = (Llam*lam*1E-10)/(h*c) # s^-1 \AA^-1

            # nlam = (Llam*912*1E-10)/(h*c) # s^-1 \AA^-1

            f = lambda l: np.interp(l, lam, nlam)
        
            n_LyC = integrate.quad(f, 10.0, 912.0)[0]
            
            n_OIII = integrate.quad(f, 10.0, lam_energy(35.11))[0]
            

            R.append(np.log10(n_OIII/n_LyC))

    
        ax.plot(grid['log10age'], np.array(R), c=color, lw=1, ls=ls)
    
    


ax.set_xlim([6.,9.])
ax.set_ylim([-4.,1.])







dummies = []
labels = []     

labels += [SPS_labels[SPS] for SPS in SPS_models]
dummies += [ax.plot([], [], ls=ls, c='k', alpha=0.5, lw = 1)[0] for ls in line_styles]        

labels += [r'$\rm Z='+str(Z)+'$' for Z in Zs]
dummies += [ax.plot([], [], ls='-', c=all_colors[all_Z==Z][0], alpha=1.0, lw = 1)[0] for Z in Zs]        


ax.legend(dummies, labels, loc = 'upper left', fontsize = 6)





ax.set_xlabel(r'$\rm\log_{10}(age/yr)$')
ax.set_ylabel(r'$\rm\log_{10}(N_{>35.11\ eV}/N_{>13.6\ eV})$')
 
    
fig.savefig('figures/hardness_OII_SPS.pdf')
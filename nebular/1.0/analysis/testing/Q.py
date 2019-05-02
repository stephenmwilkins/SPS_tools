



import numpy as np

from scipy import integrate




SPS = 'BPASSv2.binary'
IMF = 'Salpeter'

h = 6.626070040E-34 # J s

c = 3E8 # m s-1


lam = np.load('../RunGrid/inputs/'+SPS+'/lam.npy') # wavelength in \AA

nu = c/(lam*1E-10)  # frequency in Hz

Z = 0.02

for log10age in np.arange(6.0,10.0,0.2):


    input_file = SPS+'/'+'_'.join([IMF, str(Z), str(log10age)])

    # --- read in stellar SED

    Lnu = np.load('../RunGrid/inputs/'+input_file+'.npy') # (erg s^-1 Hz^-1) / M_sol yr^-1

    Llam = Lnu*nu/lam

    nu_limit = c/(91.2E-9)


    f = lambda x: np.interp(x, nu[::-1], Lnu[::-1])/(h*x)
    Lion1 = integrate.quad(f, nu_limit, nu_limit*1000)
    
    f = lambda x: np.interp(x, lam, Llam)/(h*(c*1E10/x))
    Lion2 = integrate.quad(f, 1., 912.)
    
    print log10age, np.log10(Lion1[0]) + 6., np.log10(Lion2[0]) + 6., np.log10(Lion1[0]) - np.log10(Lion2[0])


#     f = lambda x: np.interp(x, nu, Lnu)/(h*x)
#     Q1 = integrate.quad(f, nu_limit, nu_limit*100)
# 
#     f = lambda x: np.interp(x, lam, Llam)
#     Q2 = integrate.quad(f, 0., 912.)

    # print log10age, np.log10(Q1[0]) + 6., np.log10(Q2[0]) + 6.








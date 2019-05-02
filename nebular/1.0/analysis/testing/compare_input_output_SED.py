import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

c = 3E8 # m s-1

def Q(nu, Lnu):

    h = 6.626070040E-34 # J s
    nu_limit = c/(91.2E-9) # Hz (Lyman limit)
    
    f = lambda x: np.interp(x, nu[::-1], Lnu[::-1])/(h*x)
    q = integrate.quad(f, nu_limit, nu_limit*10)[0]
    
    log10Q = np.log10(q)

    return log10Q



# --------------------------------------------------------

SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_300'
Z = 0.001
age = 6.0

# --------------------------------------------------------

# --------------------------------------------------------
# --- input file


lam = np.load('../../stellar/inputs/outputs/'+SPS+'/lam.npy')
Lnu = np.load('../../stellar/inputs/outputs/'+SPS+'/'+IMF+'/'+str(Z)+'_'+str(age)+'.npy')  # W s^-1 Hz^-1

nu = c/(lam*1E-10)

print Q(nu, Lnu)

s = [(lam>921.)&(lam<1E4)]
plt.plot(np.log10(lam[s]), np.log10(Lnu[s]), lw = 3, c='0.5')

# --------------------------------------------------------
# --- output file


model = '../RunGrid/Simple/coutputs/'+SPS+'/'+IMF+'/'+str(Z)+'_'+str(age)

lam, intrinsic_stellar, stellar, nebular, linecont = np.loadtxt(model+'.cont', delimiter='\t', usecols = (0,1,2,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total

nebular -= linecont


# Cloudy input is just log10Q
# Cloudy output (for luminosity case) is initially of units of nu Lnu / (4*pi*r**2) [erg cm-2 s^-1] # C13
# Cloudy output (for luminosity case) is initially of units of nu Lnu [erg s^-1] # C17

nu = c/(lam*1E-10) # Hz

# --- change order to increasing wavelength

lam = lam[::-1]
nu = nu[::-1]

intrinsic_stellar = intrinsic_stellar[::-1]
intrinsic_stellar /= nu 
intrinsic_stellar /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
intrinsic_stellar /= 1E7 # convert to W

nebular = nebular[::-1]
nebular /= nu 
nebular /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
nebular /= 1E7 # convert to W


print Q(nu, intrinsic_stellar)




s = [(lam>921.)&(lam<1E4)]

# plt.plot(np.log10(nu[s]), np.log10(intrinsic_stellar[s]))
# plt.plot(np.log10(nu[s]), np.log10(nebular[s]))

plt.plot(np.log10(lam[s]), np.log10(intrinsic_stellar[s]))
plt.plot(np.log10(lam[s]), np.log10(nebular[s]))



plt.show()









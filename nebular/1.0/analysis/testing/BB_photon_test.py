

import numpy as np
import cPickle as pickle
from scipy import interpolate
from scipy import integrate
import read



SPS = 'BPASSv2.1.binary'
IMF = 'ModSalpeter_100'
Z = 0.001
age = 6.0

c = 3E8



model = str(Z)+'_'+str(age)

lam, intrinsic_stellar, stellar, nebular_wlines, total, linecont = np.loadtxt('../RunGrid/Simple/coutputs/'+SPS+'/'+IMF+'/'+model+'.cont', delimiter='\t', usecols = (0,1,2,3,4,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total

nu = c/(lam*1E-10) # Hz

nebular = nebular_wlines - linecont
total2 = stellar + nebular_wlines
total3 = stellar + nebular

intrinsic_stellar /= nu 
intrinsic_stellar /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
# intrinsic_stellar /= 1E7 # convert to W

total /= nu 
total /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
# total /= 1E7 # convert to W

total2 /= nu 
total2 /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
# total2 /= 1E7 # convert to W

total3 /= nu 
total3 /= 1E8 # the input was actually log10Q for a cluster log10(M)=8
# total3 /= 1E7 # convert to W



lam2 = np.arange(1000.,10000., 0.1)
nu2 = c/(lam2*1E-10) 

total4 = np.interp(nu2, nu, total3)


T = np.zeros(len(lam))
T[(lam>6400.)&(lam<6800.)] = 1.0
# T[(lam>4850)&(lam<4900)] = 1.0

T2 = np.interp(nu2, nu, T)




print integrate.trapz((1./nu) * intrinsic_stellar * T, x = nu) / integrate.trapz((1./nu) * T, x = nu), 'intrinsic stellar'

print '---------'

print integrate.trapz((1./nu) * total * T, x = nu) / integrate.trapz((1./nu) * T, x = nu), 'total'

print integrate.trapz((1./nu) * total2 * T, x = nu) / integrate.trapz((1./nu) * T, x = nu), 'stellar + nebular (inc. lines)'

L3 = integrate.trapz((1./nu) * total3 * T, x = nu) / integrate.trapz((1./nu) * T, x = nu)  # no lines

print L3, 'no lines'

print integrate.trapz((1./nu2) * total4 * T2, x = nu2) / integrate.trapz((1./nu2) * T2, x = nu2) , 'no lines, but interpolated'


lines = []
# lines += ['MgII2796','MgII2803']
# lines += ['ArIII7135','ArIII7751']
# lines += ['NeIII3869','NeIII3968','NeIV2424']
# lines += ['FeII4300','FeII2400','FeII12567','FeII16436','FeIII5271','FeIII4659','FeIV3096','FeIV2836','FeIV2829']
# lines += ['CII2325','CII2329','CII2328','CII2327','CII1335','CIII1910','CIII1907','CIV1551','CIV1548','CII1020']
# lines += ['OI6300','OI6363','OII3729','OII3726','OII2471','OII7323','OII7332','OIII1661','OIII1666','OIII5007','OIII4959','OIII2321','OVI1032','OVI1038']
# lines += ['SiII1814','SiII1308','SiII1263','SiIII1207','SiIII1892','SiIII1883','SiIV1403','SiIV1394','SiII1180']
# lines += ['NI1200','NII6584','NII6548','NII1085','NV1243','NV1239']
# lines += ['SII6720','SII4074','SII10330','SII6731','SII6716','SII4070','SIII9532','SIII9069','SIII6312','SIII3722','SIII1198','SIII1007','SIV1086']
# lines += ['HeII1640','HeII1215','HeII1085','HeII4686','HeI10830','HeI10829','HeI3889','HeI3188','HeI2945','HeI2829','HeI20581','HeI5016','HeI3965','HeI7065','HeI5876','HeI4471','HeI4026','HeI3820','HeI6678','HeI4922','HeI18685']
lines += ['HI1216','HI1026','HI6563','HI4861','HI4340','HI4102','HI3970','HI3889','HI3835','HI3798','HI18751','HI12818','HI10938','HI10049','HI9546','HI9229','HI9015','HI26252','HI21655','HI19446','HI18174']

lines = np.array(lines)


lam_lines, emergent  = read.lines(SPS, IMF, model, lines = lines)

Tint = integrate.trapz((1./nu) * T, x = nu)

for line, l,flux in zip(lines, lam_lines, emergent):

    n = 3E8/(l*1E-10)
    
    t = np.interp(n, nu, T)

    L3 += (1./n) * t * 10**flux / Tint
    
#     if t>0.0:
#     
#         print line, l, np.interp(l, lam[::-1], T[::-1]), (1./n) * t * 10**flux
    
print L3, 'lines added manually'



line_SED = np.zeros(len(lam2))


for line, l,flux in zip(lines, lam_lines, emergent):

    idx = (np.abs(lam2-l)).argmin()
    
    n = 3E8/(l*1E-10)

    line_SED[idx] += l*((10**flux)/n)/0.1

    
#     if t>0.0:
#     
#         print line, l, np.interp(l, lam[::-1], T[::-1]), (1./n) * t * 10**flux




total5 = total4 + line_SED

print integrate.trapz((1./nu2) * total5 * T2, x = nu2) / integrate.trapz((1./nu2) * T2, x = nu2) , 'lines interpolated'

print integrate.trapz((1./lam2) * total5 * T2, x = lam2) / integrate.trapz((1./lam2) * T2, x = lam2) , 'lines interpolated, lam'








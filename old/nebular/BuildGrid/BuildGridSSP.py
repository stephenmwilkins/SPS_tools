
import os
import sys
import numpy as np
import pickle
import read
from scipy import stats
from scipy import integrate

print('starting')

h = 6.626E-34
c = 3.E8

cloudy_version = 'C17.00'

version = 3.0





# --- DEFAULT
# SPS = 'BPASSv2.2.1.binary'
# IMFs = ['ModSalpeter_300']

SPSIMFs = [('BPASSv2.2.binary','ModSalpeter_300'),('P2','ModSalpeter_100')]

data_dir = '/research/astro/flare/data/SPS'




for SPSIMF in SPSIMFs:

    SPS, IMF = SPSIMF

    print('-'*50)
    print('-'*5, IMF)

    input_dir = '{data_dir}/nebular/2.0/RunGrid/Z_refQ_wdust/{cloudy_version}/{SPS}/{IMF}'.format(data_dir=data_dir, cloudy_version=cloudy_version, SPS=SPS, IMF=IMF)

    output_dir = f'{data_dir}/nebular/{version}/{SPS}/{IMF}'

    if not os.path.isdir(output_dir): os.mkdir(output_dir)

    # --- open input stellar grid (using for wavelength grid and rescaling)

    stellar_grid_filename = '{data_dir}/stellar/1.0/{SPS}/{IMF}/stellar.p'.format(data_dir=data_dir, SPS=SPS, IMF=IMF)

    stellar_grid = pickle.load(open(stellar_grid_filename,'rb'))

    log10ages = stellar_grid['log10age']

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



    SED_grid = {}

    SED_grid['Z'] = stellar_grid['Z']
    SED_grid['log10Z'] = stellar_grid['log10Z']
    SED_grid['log10age'] = log10ages
    SED_grid['lam'] = lam
    SED_grid['nu'] = 3E8/(lam*1E-10)
    SED_grid['nebular_continuum'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) ))
    SED_grid['nebular'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) ))
    # SED_grid['stellar_incident'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) )) # previously stellar
    # SED_grid['stellar_transmitted'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) ))
    SED_grid['stellar'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) )) # same as stellar_incident, included for backwards compatibility

    SED_grid['total'] = np.zeros(( len(log10ages), len(stellar_grid['Z']), len(lam) ))
    SED_grid['log10Q'] = np.zeros(( len(log10ages), len(stellar_grid['Z'])))



    line_grid = {}

    line_grid['Z'] = stellar_grid['Z']
    line_grid['log10Z'] = stellar_grid['log10Z']
    line_grid['log10age'] = log10ages
    line_grid['lines'] = lines # --- list of lines

    for line in lines:

        line_grid[line] = {}
        line_grid[line]['luminosity'] = np.zeros(( len(log10ages), len(stellar_grid['Z']) )) # no longer used
        line_grid[line]['nebular_continuum'] = np.zeros(( len(log10ages), len(stellar_grid['Z']) ))
        line_grid[line]['stellar_continuum'] = np.zeros(( len(log10ages), len(stellar_grid['Z']) ))
        # line_grid[line]['stellar_transmitted_continuum'] = np.zeros(( len(log10ages), len(stellar_grid['Z']) )) # no longer used
        line_grid[line]['total_continuum'] = np.zeros(( len(log10ages), len(stellar_grid['Z']) ))

    # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    # ---- get line fluxes and nebular contribution to continuum

    for iZ, Z in enumerate(SED_grid['Z']):

        for ia, log10age in enumerate(log10ages):



            model = '_'.join([str(Z), "{:.1f}".format(log10age)])  # --- edit this for more complicated parameters

            # ---- determine stellar log10Q

            Lnu = 1E-7 * stellar_grid['L_nu'][ia, iZ] # W s^-1 Hz^-1
            Llam = Lnu * c / (lam**2*1E-10) # W s^-1 \AA^-1
            nlam = (Llam*lam*1E-10)/(h*c) # s^-1 \AA^-1
            f = lambda l: np.interp(l, lam, nlam)
            n_LyC = integrate.quad(f, 10.0, 912.0)[0]
            log10Q_original = np.log10(n_LyC)
            SED_grid['log10Q'][ia, iZ] = log10Q_original

            # ----- continuum

            continuum_lam, continuum_nu, incident, transmitted, nebular_continuum, total, linecont = read.continuum(input_dir, model)  # --- all L_nu

            # ---- determine new log10Q

            Lnu = 1E-7 * incident # W s^-1 Hz^-1
            Llam = Lnu * c / (continuum_lam**2*1E-10) # W s^-1 \AA^-1
            nlam = (Llam*continuum_lam*1E-10)/(h*c) # s^-1 \AA^-1
            f = lambda l: np.interp(l, continuum_lam[::-1], nlam[::-1])
            n_LyC = integrate.quad(f, 10.0, 912.0)[0]
            log10Q_new = np.log10(n_LyC)


            for SED_type in ['incident', 'transmitted', 'nebular_continuum']:
                locals()[SED_type] *=10**(log10Q_original - log10Q_new)

            print(Z, log10age, log10Q_original, log10Q_new)


            # ----- lines

            line_lam, line_luminosity_intrinsic, line_luminosity_emergent  = read.lines(input_dir, model, lines = lines)

            line_luminosity_intrinsic += log10Q_original - log10Q_new
            line_luminosity_emergent += log10Q_original - log10Q_new

            line_luminosity = line_luminosity_intrinsic


            # --- build nebular line grid

            line_SED = np.zeros(len(lam))

            for l,lum in zip(line_lam, line_luminosity):

                # dl = 1.
                idx = (np.abs(lam-l)).argmin()
                dl = lam[idx] - lam[idx-1]
                n = 3E8/(l*1E-10)
                line_SED[idx] += l*((10**lum)/n)/dl

            #print(np.max(line_SED))


            SED_grid['nebular_continuum'][ia,iZ] = np.interp(lam, continuum_lam[::-1], nebular_continuum[::-1])
            SED_grid['nebular'][ia,iZ] = line_SED + SED_grid['nebular_continuum'][ia,iZ]
            SED_grid['stellar'][ia,iZ] = np.interp(lam, continuum_lam[::-1], incident[::-1]) # no dust, etc.
            SED_grid['total'][ia,iZ] = SED_grid['stellar'][ia,iZ] + SED_grid['nebular'][ia,iZ]

            for i,line in enumerate(lines):

                line_grid[line]['lam'] = line_lam[i]
                line_grid[line]['luminosity'][ia,iZ] = line_luminosity[i]
                line_grid[line]['nebular_continuum'][ia,iZ] = np.interp(line_lam[i], continuum_lam[::-1], nebular_continuum[::-1])
                line_grid[line]['stellar_continuum'][ia,iZ] = np.interp(line_lam[i], continuum_lam[::-1], incident[::-1])
                line_grid[line]['total_continuum'][ia,iZ] = line_grid[line]['nebular_continuum'][ia,iZ] + line_grid[line]['stellar_continuum'][ia,iZ]





    pickle.dump(line_grid, open(output_dir+'/lines.p', 'wb')) # --- save combined grid (can delete if necessary)
    pickle.dump(SED_grid, open(output_dir+'/nebular.p', 'wb')) # --- save combined grid (can delete if necessary)

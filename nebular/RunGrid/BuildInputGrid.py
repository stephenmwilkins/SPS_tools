



import numpy as np
import utilities as u
import os
from astropy.io import ascii
import pickle


# --- photoionisation modelling parameters

cloudy_version = 'C17.00'
log10n_H = 2.5
log10U_S_0 = -2.0
CO = 0.0
d2m = 0.3

# --- SPS/IMF choices to run

# SPS = 'BPASSv2.2.1.binary'
# IMFs = ['1p0_100','1p0_300','ModSalpeter_100', '1p7_100','1p7_300','Salpeter_100','Chabrier_100','Chabrier_300']
# IMFs = ['ModSalpeter_300']


SPS = 'BPASSv2.2.binary'
IMFs = ['ModSalpeter_300']

# SPS = 'P2'
# IMFs = ['ModSalpeter_100']




data_dir = '/research/astro/flare/data/SPS'


# --- create SPS directory if it doesn't already exist


i = 1 # model ID

for IMF in IMFs:

    output_dir = '{data_dir}/nebular/2.0/RunGrid/Z_refQ_wdust/{cloudy_version}/{SPS}/{IMF}'.format(data_dir=data_dir, cloudy_version=cloudy_version, SPS=SPS, IMF=IMF)

    print(output_dir)

    if not os.path.exists(output_dir): os.mkdir(output_dir)


    # --- read in stellar grid

    stellar_grid_filename = '{data_dir}/stellar/1.0/{SPS}/{IMF}/stellar.p'.format(data_dir=data_dir, SPS=SPS, IMF=IMF)

    stellar_grid = pickle.load(open(stellar_grid_filename,'rb'))

    print(stellar_grid['Z'])
    print(stellar_grid['log10Z'])
    print(stellar_grid['log10age'])
    print(stellar_grid['L_nu'].shape)

    lam = stellar_grid['lam']

    REF_iZ = (np.abs(stellar_grid['log10Z'] - (-2.0))).argmin()

    log10Q_REF = u.measure_log10Q(lam, stellar_grid['L_nu'][0, REF_iZ])

    summary = {'Z': [], 'log10age':[], 'log10U_S': [], 'log10Q': [], 'log10Q_orig': []}

    for iZ, Z in enumerate(stellar_grid['Z']):

        for ia, log10age in enumerate(stellar_grid['log10age']):

            CLOUDY_SED = u.create_CLOUDY_SED(lam, stellar_grid['L_nu'][ia, iZ])

            log10Q_orig = u.measure_log10Q(lam, stellar_grid['L_nu'][ia, iZ])

            log10U_S = log10U_S_0 + (log10Q_orig - log10Q_REF)/3.

            log10Q = u.determine_log10Q(log10U_S, log10n_H) # --- this is the input

            print('{0} {1:.2f} | {2:.2f} {3:.2f} {4:.2f}'.format(Z, log10age, log10U_S, log10Q, log10Q_orig))

            summary['Z'].append(Z)
            summary['log10age'].append(log10age)
            summary['log10U_S'].append(log10U_S)
            summary['log10Q'].append(log10Q)
            summary['log10Q_orig'].append(log10Q_orig)

            # ----- start CLOUDY input list

            cinput = []

            # --- Define the incident radiation field shape

            cinput.append('interpolate{      10.0000     -99.0000}\n')
            for i1 in range(int(len(CLOUDY_SED)/5)): cinput.append('continue'+''.join(CLOUDY_SED[i1*5:(i1+1)*5])+'\n')

            # --- Define the chemical composition

            a = u.abundances(Z, CO, d2m) # --- determine abundances

            a_nodep =  u.abundances(Z, CO, 0.0) # --- determine abundances for no depletion

            for ele in ['He'] + u.metals:

                cinput.append('element abundance '+u.name[ele]+' '+str(a[ele])+'\n')

            # cinput.append('element off limit -7') # should speed up the code


            # graphite, scale by total C abundance relatve to ISM
            scale = 10**a_nodep['C']/2.784000e-04
            cinput.append(f'grains Orion graphite {scale}'+'\n')

            # silicate, scale by total Si abundance relatve to ISM NOTE: silicates also rely on other elements.
            scale = 10**a_nodep['Si']/3.280000e-05
            cinput.append(f'grains Orion silicate {scale}'+'\n')


            # --- Define the luminosity

            cinput.append('Q(H) '+str(log10Q)+'\n')

            # --- Define the geometry

            cinput.append('hden '+str(log10n_H)+' log constant density\n')
            cinput.append('radius -2 parsecs\n')
            cinput.append('sphere\n')
            cinput.append('covering factor 1.0 linear\n')

            # --- Processing commands

            cinput.append('iterate to convergence\n')
            cinput.append('set temperature floor 100 linear\n')
            cinput.append('stop temperature 4000K\n')
            cinput.append('stop efrac -2\n')


            # --- define output filename


            output_file = '_'.join([str(Z), "{:.1f}".format(log10age)])

            cinput.append('save last continuum "'+output_dir+'/'+output_file+'.cont" units Angstroms no clobber\n')
            cinput.append('save last lines, array "'+output_dir+'/'+output_file+'.lines" units Angstroms no clobber\n')
            cinput.append('save overview "'+output_dir+'/'+output_file+'.ovr" last\n')

            # --- write cloudy input file

            f = open('cinputs/'+str(i)+'.in','w')
            f.writelines(cinput)
            f.close()


            i += 1


    ascii.write(summary,'input_summary_'+IMF+'.dat', names=['Z', 'log10age', 'log10U_S', 'log10Q', 'log10Q_orig'])
    print('created '+str(i-1)+' models')

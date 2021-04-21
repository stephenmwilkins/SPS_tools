

import numpy as np
from scipy import integrate
import math
import pickle
import h5py
import yaml

from . import abundances


parameter_names = ['SPS', 'IMF', 'ia', 'iZ', 'cloudy_version', 'log10n_H', 'log10U_S_0', 'CO', 'd2m', 'REF_log10Z', 'REF_log10age', 'log10radius', 'covering_factor', 'stop_T', 'stop_efrac', 'T_floor', 'stellar_grid_file', 'output_file']



def default_params():

    # ---------------------------------------------------
    # --- SPS SED parameters

    SPS = 'BPASSv2.2.1.binary'
    IMF = 'ModSalpeter_300'
    ia = 0 # age index of the SED to use, range/mapping to Myr depends on choice of SPS model
    iZ = 0 # metallicity index of the SED to use, range/mapping to Z depends on choice of SPS model

    # ---------------------------------------------------
    # --- Photoionisation parameters

    cloudy_version = 'C17.00' # cloudy version

    # --- core parameters

    log10n_H = 2.5 # Hydrogen density
    log10U_S_0 = -2.0 # ionisation parameter at the reference metallicity and age ## NOTE: this is not how other people handle this but I think it makes more sense
    CO = 0.0 #
    d2m = 0.3 # dust-to-metal ratio

    # --- other parameters

    REF_log10Z = -2.0 # reference metallicity
    REF_log10age = 6.0 # reference age
    log10radius = -2.0 # radius in log10 parsecs
    covering_factor = 1.0 # covering factor. Keep as 1 as it is more efficient to simply combine SEDs to get != 1.0 values

    # --- Processing commands

    stop_T = 4000 # K
    stop_efrac = -2
    T_floor = 100 # K


    # ---------------------------------------------------
    # --- file locations

    # --- inputs
    stellar_grid_file = None

    # --- outputs
    output_file = 'test'


    # ---------------------------------------------------
    # --- return a dictionary of these parameters

    params = {}
    for name in parameter_names:
        params[name] = locals()[name]

    return params



def build_input(SPS, IMF, ia, iZ, cloudy_version, log10n_H, log10U_S_0, CO, d2m, REF_log10Z, REF_log10age, log10radius, covering_factor, stop_T, stop_efrac, T_floor, stellar_grid_file, output_file):

    # --- this function builds the cloudy input file



    # --- read in pure stellar grid

    if stellar_grid_file.split('.')[-1] == 'p':
        stellar_grid = pickle.load(open(stellar_grid_file,'rb'))
        lam = stellar_grid['lam']

    if stellar_grid_file.split('.')[-1] == 'h5':
        # stellar_grid = pickle.load(open(stellar_grid_file,'rb'))
        lam = stellar_grid['lam']


    # --- determine the metallicity and age indicies for the reference metallicity and age

    iZ_REF = (np.abs(stellar_grid['log10Z'] - (REF_log10Z))).argmin()
    ia_REF = (np.abs(stellar_grid['log10age'] - (REF_log10age))).argmin()

    # --- calculate the LyC flux for the reference metallicity and age SED
    log10Q_REF = measure_log10Q(lam, stellar_grid['L_nu'][ia_REF, iZ_REF])

    # --- calculate the ionising photon luminosity for the target SED
    log10Q_orig = measure_log10Q(lam, stellar_grid['L_nu'][ia, iZ])

    # --- calculate the actual ionisation parameter for the target SED. Only when the target SED == reference SED will log10U_S = log10U_S_0
    log10U_S = log10U_S_0 + (log10Q_orig - log10Q_REF)/3.

    # --- now determine the actual ionising photon luminosity for the target SED. This is the normalising input to CLOUDY
    log10Q = determine_log10Q(log10U_S, log10n_H)

    # --- convert the SED to the CLOUDY format. The normalisation doesn't matter as that is handled by log10Q above
    CLOUDY_SED = create_CLOUDY_SED(lam, stellar_grid['L_nu'][ia, iZ])

    # --- get metallicity
    Z = stellar_grid['Z'][iZ]

    # --- determine elemental abundances for given Z, CO, d2m, with depletion taken into account
    a = abundances.abundances(Z, CO, d2m)

    # --- determine elemental abundances for given Z, CO, d2m, WITHOUT depletion taken into account
    a_nodep =  abundances.abundances(Z, CO, 0.0) # --- determine abundances for no depletion


    # ----- start CLOUDY input file (as a list)
    cinput = []

    # --- Define the incident radiation field shape
    cinput.append('interpolate{      10.0000     -99.0000}\n')
    for i1 in range(int(len(CLOUDY_SED)/5)): cinput.append('continue'+''.join(CLOUDY_SED[i1*5:(i1+1)*5])+'\n')

    # --- Define the chemical composition
    for ele in ['He'] + abundances.metals:
        cinput.append('element abundance '+abundances.name[ele]+' '+str(a[ele])+'\n')

    # cinput.append('element off limit -7') # should speed up the code

    # --- add graphite and silicate grains
    # graphite, scale by total C abundance relatve to ISM
    scale = 10**a_nodep['C']/2.784000e-04
    cinput.append(f'grains Orion graphite {scale}'+'\n')
    # silicate, scale by total Si abundance relatve to ISM NOTE: silicates also rely on other elements.
    scale = 10**a_nodep['Si']/3.280000e-05
    cinput.append(f'grains Orion silicate {scale}'+'\n')


    # --- Define the luminosity
    cinput.append(f'Q(H) {log10Q}\n')

    # --- Define the geometry

    cinput.append(f'hden '+str(log10n_H)+' log constant density\n')
    cinput.append(f'radius {log10radius} parsecs\n')
    cinput.append(f'sphere\n')
    cinput.append(f'covering factor {covering_factor} linear\n')

    # --- Processing commands

    cinput.append(f'iterate to convergence\n')
    cinput.append(f'set temperature floor {T_floor} linear\n')
    cinput.append(f'stop temperature {stop_T}K\n')
    cinput.append(f'stop efrac {stop_efrac}\n')

    # --- define output filename

    cinput.append(f'save last continuum "{output_file}.cont" units Angstroms no clobber\n')
    cinput.append(f'save last lines, array "{output_file}.lines" units Angstroms no clobber\n')
    cinput.append(f'save overview "{output_file}.ovr" last\n')

    # --- write input file
    open(f'{output_file}.in','w').writelines(cinput)


    # --- make a dictionary of the parameters including derived parameters

    derived_parameter_names = ['log10Q_orig', 'log10U_S', 'log10Q', 'Z']

    params = {}
    for name in parameter_names + derived_parameter_names:
        params[name] = locals()[name]

    # --- write parameter file, including derived parameters
    yaml.dump(params, open(f'{output_file}.yaml','w'))

    return cinput




def determine_log10Q(log10U_S, log10n_H):


    alpha_B = 2.59E-13 # cm3 s-1
    c = 3.0E8 # m s-1
    c_cm = c * 100. # cm s-1

    n_H = 10**log10n_H
    U_S = 10**log10U_S

    epsilon = 1.

    Q = ((U_S*3.*c_cm)**3/alpha_B**2)*((4.*np.pi)/(3*epsilon**2*n_H))

    return np.log10(Q)




def measure_log10Q(lam, Lnu):

    h = 6.626070040E-34 # J s
    c = 3E8 # m s-1
    nu = c/(lam * 1E-10) # nu in Hz assuming lam in \AA
    nu_limit = c/(91.2E-9)
    f = lambda x: np.interp(x, nu[::-1], Lnu[::-1])/(h*x)
    log10Q = np.log10(integrate.quad(f, nu_limit, nu_limit*1000)[0]) # 1 M_sol mass star cluster

    return log10Q


def create_CLOUDY_SED(lam, Lnu):


    nu = 3E8/(lam*1E-10)  # frequency in Hz
    nu_log10 = np.log10(nu)

    Lnu_log10 = np.log10(Lnu+1E-99)
    Lnu_log10 -= np.max(Lnu_log10)

    # --- reformat for CLOUDY

    CLOUDY_SED = ['{'+'{0:.5f} {1:.5f}'.format(x,y)+'}' for x,y in zip(nu_log10[::2], Lnu_log10[::2])]
    CLOUDY_SED = CLOUDY_SED[:19000]
    CLOUDY_SED = CLOUDY_SED[::-1]

    return CLOUDY_SED

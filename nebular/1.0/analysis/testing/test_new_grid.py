

import cPickle as pickle
import numpy as np
from scipy import ndimage



def read_grid(data_dir, grid, lines):
    
    R = {}
    R['lam'] = {}
    R['fluxes'] = {}
    R['nebular_continuum'] = {}
    R['stellar_continuum'] = {}
    
    first_instance = True
    
    for line in lines:
        for l in line.split(','):
        
            nR = pickle.load(open(data_dir + grid+'/LINES/'+l+'.p', 'r'))
                
            R['lam'][l] = nR['lam'] 
            R['fluxes'][l] = nR['line_flux']
            R['nebular_continuum'][l] = nR['nebular_continuum']
            R['stellar_continuum'][l] = nR['stellar_continuum'] 
               
            if first_instance:               
                for p in ['log10age', 'Z', 'log10n_H', 'log10U', 'd2m', 'CO']: R[p] = nR[p]
                R['log10Z'] = np.log10(R['Z'])                 
                first_instance = False

    return R
    



def getEW(R, p, line):

    interp_order = 1

    nebular_params = [[np.interp(p[parameter], R[parameter], range(len(R[parameter])))] for parameter in ['Z','log10n_H','log10U','d2m','CO']]
    stellar_params = [[np.interp(p[parameter], R[parameter], range(len(R[parameter])))] for parameter in ['log10age','Z']]

    lam = np.mean([R['lam'][l] for l in line.split(',')])
    flux = np.sum([10**ndimage.map_coordinates(R['fluxes'][l], nebular_params, order=interp_order)[0] for l in line.split(',')])
    nebular_continuum = np.mean([10**ndimage.map_coordinates(R['nebular_continuum'][l], nebular_params, order=interp_order)[0] for l in line.split(',')])
    stellar_continuum = np.mean([10**ndimage.map_coordinates(R['stellar_continuum'][l], stellar_params, order=interp_order)[0] for l in line.split(',')])
    
    EW = np.log10(R['lam'][l]*flux/(nebular_continuum+stellar_continuum))

    print flux, nebular_continuum, stellar_continuum, '|', 10**EW

    return EW






lines = ['CIII1907,CIII1910']

p = {}
p['Z'] = 0.004
p['log10age'] = 7.
p['CO'] = 0.0
p['d2m'] = 0.3
p['log10n_H'] = 2.0
p['log10U'] = -2.5
 

grid = 'BPASSv2.binary_Salpeter'

# ---- new grid

data_dir = '../BuildGrid/grids/'
R = read_grid(data_dir, grid, lines)

getEW(R, p, lines[0])



# ---- old grid

data_dir = '../../../PhotoionisationGrid_OLD/0.10/BuildGridNew/grids/'
R = read_grid(data_dir, grid, lines)

getEW(R, p, lines[0])





import numpy as np







def lines(SPS, IMF, model, lines = [], grid = 'Simple', cloudy_version = 'C17.00'):

    data_dir = '../RunGrid/'+grid+'/coutputs/'+cloudy_version+'/'
    
    log10M_cluster = 8.


    lam, cID, intrinsic, emergent = np.loadtxt(data_dir + SPS+ '/'+IMF+'/'+model+'.lines', dtype = str, delimiter='\t', usecols = (0,1,2,3)).T

    lam = lam.astype(float) 
    intrinsic = intrinsic.astype(float)        
    emergent = emergent.astype(float) 

    Lam = []
    ID = []
    Emergent = []


    for l,id, intrin, emerg in zip(lam, cID, intrinsic, emergent): 

        l = int(np.round(l,0))

        li = filter(None, id.split(' '))  

        e = li[0]

        i = li[1]
        j = '-'
        if i == '1': j = 'I'
        if i == '2': j = 'II'
        if i == 'II': j = 'II'
        if i == '3': j = 'III'
        if i == '4': j = 'IV'
        if i == '5': j = 'V'
        if i == '6': j = 'VI'
        if i == '7': j = 'VII'
        if i == '8': j = 'VIII'
        if i == '9': j = 'IX'
        if i == '10': j = 'X'
        if i == '11': j = 'XI'
        if i == '12': j = 'XII'
        if i == '13': j = 'XIII'
        if i == '14': j = 'XIV'


        nid = e+j+str(l)
    
        Lam.append(l)
        ID.append(nid)
        Emergent.append(emerg)
            
                
    Lam = np.array(Lam)
    ID = np.array(ID)
    Emergent = np.array(Emergent)

    n_Lam = []
    n_emergent = []

    for line in lines:
    
        if line in ID:
        
            n_emergent.append(Emergent[(ID==line)][0])
            n_Lam.append(Lam[(ID==line)][0])
            
        else:
        
            n_emergent.append(-99.)
            n_Lam.append(1000.)


    n_Lam = np.array(n_Lam)
    n_emergent = np.array(n_emergent) - log10M_cluster # correct for size of cluster # erg s^-1

    return n_Lam, n_emergent







def nebular_continuum(SPS, IMF, model, grid = 'Simple', cloudy_version = 'C17.00'):

    data_dir = '../RunGrid/'+grid+'/coutputs/'+cloudy_version+'/'
    
    log10M_cluster = 8.

    # ----- Open SED
    
    continuum_lam, continuum_wlines, linecont = np.loadtxt(data_dir + SPS+'/'+IMF+'/'+model+'.cont', delimiter='\t', usecols = (0,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
  
    continuum = continuum_wlines - linecont # remove contribution of lines (necessary for equivalent width measurements

#     continuum_wlines = 4.*np.pi * (100.*0.01*3.0857E16)**2 # L erg s^-1 # no longer needed in C17
#     continuum = 4.*np.pi * (100.*0.01*3.0857E16)**2 # L erg s^-1 # no longer needed in C17

    continuum_wlines /= 10**log10M_cluster
    continuum /= 10**log10M_cluster


    return continuum_lam, continuum, continuum_wlines






def stellar_continuum(SPS, IMF, model, grid = 'Simple', cloudy_version = 'C17.00'):

    data_dir = '../RunGrid/'+grid+'/coutputs/'+cloudy_version+'/'
    
    log10M_cluster = 8.

    # ----- Open SED
    
    continuum_lam, stellar = np.loadtxt(data_dir + SPS+'/'+IMF+'/'+model+'.cont', delimiter='\t', usecols = (0,1)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
  
    stellar /= 10**log10M_cluster

    return continuum_lam, stellar

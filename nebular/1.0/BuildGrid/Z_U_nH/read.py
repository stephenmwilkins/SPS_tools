
import numpy as np







def lines(SPS, IMF, model, lines = [], cloudy_grid = 'Z', cloudy_version = 'C17.00'):

    data_dir = '../../RunGrid/'+cloudy_grid+'/coutputs/'+cloudy_version+'/'
    
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

        li = list(filter(None, id.split(' '))) 

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




def continuum(SPS, IMF, model, cloudy_grid = 'Z', cloudy_version = 'C17.00'):

    data_dir = '../../RunGrid/'+cloudy_grid+'/coutputs/'+cloudy_version+'/'
    
    
    log10M_cluster = 8.

    # ----- Open SED
    
    lam, stellar, nebular, total, linecont  = np.loadtxt(data_dir + SPS+'/'+IMF+'/'+model+'.cont', delimiter='\t', usecols = (0,1,3,4,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total
  
  
    nu = 3E8/(lam*1E-10)
  
    nebular_continuum = nebular - linecont
  
  
    for SED_type in ['stellar', 'nebular_continuum', 'total', 'linecont']:
  
        locals()[SED_type] /= 10**log10M_cluster
        locals()[SED_type] /= nu
        


    return lam, nu, stellar, nebular_continuum, total, linecont
    
    
    
    
    
    
    

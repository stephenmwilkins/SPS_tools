
import numpy as np

def nebular(SPS, model, lines = []):

    data_dir = '../RunGrid/coutputs/'
    
    log10M_cluster = 6.

    # ----- Open SED
    
    SEDlam, stellar, nebular, linecont = np.loadtxt(data_dir + SPS+ '/'+model+'.cont', delimiter='\t', usecols = (0,2,3,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total

    nebular -= linecont # remove contribution of lines (necessary for equivalent width measurements
    
    stellar *= 4.*np.pi * (100.*0.01*3.0857E16)**2 # L erg s^-1 
    nebular *= 4.*np.pi * (100.*0.01*3.0857E16)**2 # L erg s^-1 

    log10SEDlam = np.log10(SEDlam)

    nebular = np.log10(nebular) - log10M_cluster # correct for size of cluster
    stellar = np.log10(stellar) - log10M_cluster # correct for size of cluster



    # --------------------------------------------------------
    # ----- Lines

    lam, cID, intrinsic, emergent = np.loadtxt(data_dir + SPS+ '/'+model+'.lines', dtype = str, delimiter='\t', usecols = (0,1,2,3)).T

    lam = lam.astype(float) 
    loglam = np.log10(lam)
    intrinsic = intrinsic.astype(float)        
    emergent = emergent.astype(float) 

    Lam = []
    element = []
    ion = []
    ID = []
    Intrinsic = []
    Emergent = []

    continuum_nebular = []
    continuum_stellar = []

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
    
        if nid in lines:
        
            Lam.append(l)
            ID.append(nid)
            element.append(e)
            ion.append(j)
            Intrinsic.append(intrin)
            Emergent.append(emerg)
            
            continuum_nebular.append( np.interp(l, SEDlam[::-1], nebular[::-1]) )
            continuum_stellar.append( np.interp(l, SEDlam[::-1], stellar[::-1]) )

    Lam = np.array(Lam)
    element = np.array(element)
    ion = np.array(ion)
    ID = np.array(ID)

    Intrinsic = np.array(Intrinsic) - log10M_cluster # correct for size of cluster # erg s^-1
    Emergent = np.array(Emergent) - log10M_cluster # correct for size of cluster # erg s^-1

    continuum_nebular = np.array(continuum_nebular)
    continuum_stellar = np.array(continuum_stellar)

    return Lam, ID, Emergent, continuum_nebular, continuum_stellar



    
    
def stellar(SPS, model, line_IDs, line_lams):

    c = 3E8

    data_dir = '../RunGrid/inputs/'

    lam = np.load(data_dir + SPS+'/lam.npy')

    Lnu = np.load(data_dir + SPS+'/' + model + '.npy') # W s^-1 Hz^-1
    
    nu = c/(lam*1E-10)
    
    L = Lnu*1E7*nu # erg s^-1
      
    stellar_continuum = np.array([np.interp(l, lam, L) for l in line_lams])
    
    return stellar_continuum
    
        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
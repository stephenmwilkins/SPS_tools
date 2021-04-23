
import numpy as np







def lines(output_dir, output_file, lines = []):


    log10M_cluster = 8.

    lam, cID, intrinsic, emergent = np.loadtxt(f'{output_dir}/{output_file}.lines', dtype = str, delimiter='\t', usecols = (0,1,2,3)).T

    lam = lam.astype(float)
    intrinsic = intrinsic.astype(float)
    emergent = emergent.astype(float)

    Lam = []
    ID = []
    Intrinsic = []
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
        Intrinsic.append(intrin)


    Lam = np.array(Lam)
    ID = np.array(ID)
    Emergent = np.array(Emergent)
    Intrinsic = np.array(Intrinsic)

    n_Lam = []
    n_emergent = []
    n_intrinsic = []

    for line in lines:

        if line in ID:

            n_emergent.append(Emergent[(ID==line)][0])
            n_intrinsic.append(Intrinsic[(ID==line)][0])
            n_Lam.append(Lam[(ID==line)][0])

        else:

            n_emergent.append(-99.)
            n_intrinsic.append(-99.)
            n_Lam.append(1000.)


    n_Lam = np.array(n_Lam)
    n_emergent = np.array(n_emergent) - log10M_cluster # correct for size of cluster # erg s^-1
    n_intrinsic = np.array(n_intrinsic) - log10M_cluster # correct for size of cluster # erg s^-1


    return n_Lam, n_intrinsic, n_emergent




def continuum(output_dir, output_file):

    # ----- Open SED

    lam, incident, transmitted, nebular, total, linecont  = np.loadtxt(f'{output_dir}/{output_file}.cont', delimiter='\t', usecols = (0,1,2,3,4,8)).T # 1 = incident, 2 = transmitted, 3 = nebular, 4 = total, 8 = contribution of lines to total

    # --- frequency
    nu = 3E8/(lam*1E-10)

    # --- nebular continuum is the total nebular emission (nebular) - the line continuum (linecont)
    nebular_continuum = nebular - linecont


    for SED_type in ['incident', 'transmitted', 'nebular_continuum', 'total', 'linecont']:

        locals()[SED_type] /= 10**7  # 
        locals()[SED_type] /= nu



    return lam, nu, incident, transmitted, nebular_continuum, total, linecont

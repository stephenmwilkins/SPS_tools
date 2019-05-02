


import numpy as np
import read
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm



plt.rcParams['mathtext.fontset'] = 'stixsans'
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 8

plt.rcParams['ytick.labelsize'] = 6 
plt.rcParams['xtick.labelsize'] = 6 
 
plt.rcParams['ytick.direction'] = 'in'    # direction: in, out, or inout
plt.rcParams['xtick.direction'] = 'in'    # direction: in, out, or inout

plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.minor.visible'] = True    



fig = plt.figure(figsize=(5,3 ))

left = 0.1
width = 0.8
bottom = 0.15
height = 0.75
legend = 0.1

ax = fig.add_axes([left, bottom, width, height])

ax_legend = fig.add_axes([left+width, bottom, legend, height])

lam_limits = [3., 3.55]
R_limit =  -2.0 



ax.axhline(0.0, lw = 2, alpha = 0.1, c= 'k')


# -------------------

default = True

if default:

    SPS = 'BPASSv2.1.binary'
    IMF = 'ModSalpeter_300'
    
    SPSIMF = 'BPASSv2.1\quad Modified\, Salpeter\,(m=0.1-300\,M_{\odot})'
    
    age = 6.0
    Z = 0.004



model = '_'.join([str(Z), str(age)]) 




lines = read.BASIClines(SPS, IMF, model, grid = 'Simple')

sel = [(np.log10(lines.lam)>lam_limits[0])&(np.log10(lines.lam)<lam_limits[1])&(lines.R>R_limit)]

lam = lines.lam[sel]
ID = lines.ID[sel]
element = lines.element[sel]
R = lines.R[sel]


dc = cm.gray(0.5)
c = {}
c['O'] = cm.rainbow(0.2)
c['H'] = cm.rainbow(1.0)
c['N'] = cm.rainbow(0.4)
c['C'] = cm.rainbow(0.0)
c['S'] = cm.rainbow(0.6)
c['He'] = cm.rainbow(0.8)
c['Si'] = cm.viridis(0.3)
c['Fe'] = dc
c['Ti'] = dc
c['Ne'] = dc
c['Ar'] = dc
c['Na'] = dc
c['Mg'] = dc
c['Li'] = dc
c['B'] = dc
c['Cr'] = dc
c['Cl'] = dc
c['Ca'] = dc
c['blend'] = cm.gray(0.3)




include_elements = ['H','He','S','N','Si','O','C','Fe', 'Ti', 'Ne', 'Ar', 'Na', 'Mg', 'Li', 'B', 'Cr', 'Cl', 'blend'][::-1]
ie = 0

for e in include_elements:

    s = [(element==e)]

    

    if len(ID[s])>0:

        print "lines += ['"+"','".join(ID[s])+"']"

        if e == 'blend':
            ax.scatter(lam[s], R[s], c = [c[e]]*len(ID[s]), lw = 0.3, alpha = 0.5, s = 20, hatch = 5*'//')
        else: 
            ax.scatter(lam[s], R[s], c = [c[e]]*len(ID[s]), lw = 0, alpha = 0.5, s = 20)



        for i, id in enumerate(ID[s]):
            ax.annotate(id, (lam[s][i]+25, R[s][i]-0.037), size = 4, alpha = 0.6, family = 'sans-serif')

            # print id, element[s][i]

        
        
        # --- add to legend
        
        if e == 'blend':
            ax_legend.scatter(0.2, ie, c = c[e], lw = 0.3, alpha = 0.5, s = 20, hatch = 5*'//')
        else: 
            ax_legend.scatter(0.2, ie, c = c[e], lw = 0, alpha = 0.5, s = 20)
        
        ax_legend.text(0.3, ie - 0.05, e, color = '#555555', weight = 'medium', fontsize = 5, ha='left', va='center')

        ie += 1
    

# 
# # --- read in stack list
# 
# include_obs = True
# 
# if include_obs:
# 
#     stack = np.loadtxt('stack_linelist.dat', dtype = str).T
# 
#     for id in stack[2]:
# 
#         try:
#             s = [ID==id]
# 
#             l = np.log10(lam[s])[0]
#             r = R[s][0]
#             e = element[s][0]
# 
#             print id, l,r,e
#     
#             ax.scatter(l, r, c = np.array([c[e]]), lw = 0.5, edgecolor = 'k', s = 20)
#         except:
#     
#             print id, 'didnt find line'
# 
# 



# --- make legend

ax_legend.set_xlim([0.,1])
ax_legend.set_ylim([ax_legend.set_ylim()[0]-6, ax_legend.set_ylim()[1]+6])
ax_legend.set_axis_off()



ax.minorticks_on()


# ax.set_xlabel(r"${\rm\lambda/\AA}$")
# ax.set_ylabel(r"${\rm\log_{10}(L_{i}/L_{H\beta})}$")

ax.set_xlabel(r"${\rm Wavelength/\AA}$")
ax.set_ylabel(r"${\rm\log_{10}(line\ luminosity\ relative\ to\ H\beta)}$")

# ax.set_xlabel(r"${\rm\log_{10}(\lambda/\AA)}$")




model_label = '\quad'.join([SPSIMF, 'Cloudy\,Version=C17.00', 'Z='+str(Z), r'\log_{10}(n/cm^{-3})=2.0', r'\log_{10}(U_S)=-2.5', r'\log_{10}[(C/O)/(C/O)_{\odot}]=0.0', r'\xi_{d}=0.3'])
plt.figtext(left+width/2., bottom+height+0.01, r'${\rm '+model_label+'}$', size = 5, alpha = 0.3, ha = 'center', va = 'bottom')





if default:
    fig_name = 'figures/UVlines_'+str(Z)+'.pdf'
else:
    fig_name = 'figures/UVlines_'+model+'.pdf'
    
    
print fig_name

fig.savefig(fig_name)
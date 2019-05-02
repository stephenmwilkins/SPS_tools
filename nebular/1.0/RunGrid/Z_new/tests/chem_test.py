import numpy as np
import utilities as u



Z = 0.0131 # solar

CO = 0.0

d2m = 0.0





print '------------'

print 'Z_sol = ', u.metallicity(u.sol)

print '------------'

# determine solar dust depletion factor

print 'solar dust depletion factor:', u.dust_to_metal(u.sol, u.depsol)


print '------------'

a = u.abundances(Z, CO, d2m)

for e in a.keys(): print e, a[e], u.sol[e], u.depsol[e]

print '------------ CO'

print a['C'] - a['O']

print '------------ CO'

print u.mf(Z)

print u.metallicity(a)
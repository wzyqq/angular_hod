import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmos = FlatLambdaCDM(H0=100, Om0=0.3)
z_high = 1.2
z_low  = 1.1
distance_high = cosmos.comoving_distance(z_high).value
distance_low = cosmos.comoving_distance(z_low).value
volumn = 4.0/3*np.pi*(distance_high**3-distance_low**3)

observe = len(np.loadtxt("../../xuhai_cat/cosmos_ngdata"))
ng_average = observe/(volumn*2.06/41252.96125)
ng_std = ng_average/5.0
print(ng_average,ng_std)
#np.savetxt("ng_average",ng_average,fmt='%14.8lf')
#np.savetxt("ng_std",ng_std,fmt='%14.8lf')

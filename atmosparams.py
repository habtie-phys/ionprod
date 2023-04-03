import numpy as np
from pymsis import msis
import pymsis.utils as utlis
def get_params(dt, alts_km, glat, glon):
    (f107, f107a, ap) = utlis.get_f107_ap(dt)
    aps = [list(ap)]
    atmos =  msis.run(dt, glon, glat, alts_km, f107, f107a, aps, geomagnetic_activity = -1, version=2.1)
    atmos = np.squeeze(atmos)
    #
    rho = atmos[:, 0]
    nN2 = atmos[:, 1]
    nO2 = atmos[:, 2]
    nO = atmos[:, 3]
    nHe = atmos[:, 4]
    nH = atmos[:, 5]
    nAr = atmos[:, 6]
    nN = atmos[:, 7]
    nNO = atmos[:, 9]
    Tn = atmos[:, 10]
    #
    knan = np.isnan(nH)
    nH[knan] = 1.
    knan = np.isnan(nHe)
    nHe[knan] = 1.
    knan = np.isnan(nN)
    nN[knan] = 1.
    knan = np.isnan(nNO)
    nNO[knan] = 1.
    knan = np.isnan(nO)
    nO[knan] = 1.
    #
    z = alts_km*1e3 # km to m
# Physical constants
    kB	= 1.380662e-23;		    # Boltzmann constant [J/K]
    Re	= 6.378e6;		    # Radius of earth [m]

    # Altitude varying gravitational acceleration
    g = 9.82*(Re/(Re+z))**2
    #Average molecular mass (kg) and scale height (cm)
    m = rho/(nO + nO2 + nN2 + nN + nHe + nH + nAr)
    H = (kB*Tn/(m*g))*1e2
    rho_g_cm3 = 1e-3*rho # in g/cm^3
    #
    return (rho_g_cm3, m, H)
# %%

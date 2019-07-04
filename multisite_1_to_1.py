import numpy as np
from scipy.integrate import odeint
def multisite_1_to_1(a, x, kdax):
    def ode_multisite_1_to_1(concs, t, kdax):
        KON = 1
        KOFF = kdax/KON
        p = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        l = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        pl = KON*concs[0]*concs[1]-(KOFF)*concs[2]
        return [p, l, pl]

    return odeint(ode_multisite_1_to_1, [a, x, 0], np.linspace(0, 1, 100), args=(kdax,))[-1, 2]/x

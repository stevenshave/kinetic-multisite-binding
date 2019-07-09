import numpy as np
from scipy.integrate import solve_ivp, odeint
from multisite_1_to_1 import multisite_1_to_1


import numpy as np
import matplotlib.pyplot as plt

def multisite_1_to_1_solve_ivp(a, x, kd1, interval=(0,1)):
    def ode_multisite_1_to_1(t, concs, kd1):
        a, x, a1_x = concs
        r1 = -a*x + kd1*a1_x
        dadt = 0.0 +r1 
        da1_xdt = 0.0 -r1 
        dxdt = r1
        return [dadt, dxdt, da1_xdt]
    res = solve_ivp(lambda t, y: ode_multisite_1_to_1(t,y, kd1), interval, [a, x,0.0], method="LSODA").y[2,-1]/x
    return res





import numpy as np
from scipy.integrate import odeint
def multisite_1_to_1_odeint(a, x, kdax):
    def ode_multisite_1_to_1(concs, t, kdax):
        KON = 1
        KOFF = kdax/KON
        p = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        l = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        pl = KON*concs[0]*concs[1]-(KOFF)*concs[2]
        return [p, l, pl]

    return odeint(ode_multisite_1_to_1, [a, x, 0.0],[0,1,2,3,4,5,6,7,8,9,10], args=(kdax,))[-1,2]/x

print(multisite_1_to_1_solve_ivp(20.0, 20.0, 10.0))
print(multisite_1_to_1_odeint(20.0, 20.0, 10.0))

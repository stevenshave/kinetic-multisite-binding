import numpy as np
from scipy.integrate import solve_ivp
def multisite_1_to_2(a, x, kd1, kd2, interval=(0,1)):
	def ode_multisite_1_to_2(concs, t, kd1, kd2):
		a, x, a1_x, a2_x, a1_2_x = concs
		r1 = -a*x + kd1*a1_x
		r2 = -a*x + kd2*a2_x
		r3 = -a1_x*x + kd2*a1_2_x
		r4 = -a2_x*x + kd1*a1_2_x
		dadt = 0.0 +r1 +r2 
		da1_xdt = 0.0 -r1 +r3 
		da2_xdt = 0.0 -r2 +r4 
		da1_2_xdt = 0.0 -r3 -r4 
		dxdt = r1+r2+r3+r4
		return [dadt, dxdt, da1_xdt, da2_xdt, da1_2_xdt]
	res = solve_ivp(lambda t, y: ode_multisite_1_to_2(y,t,kd1, kd2), interval, [a, x,0.0, 0.0, 0.0]).y[2:,-1]
	return (1*(sum(res[0:2]))+2*res[2])/x

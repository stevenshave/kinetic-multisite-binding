import numpy as np
from scipy.integrate import solve_ivp
def multisite_1_to_1(a, x, kd1, interval=(0,1)):
	def ode_multisite_1_to_1(concs, t, kd1):
		a, x, a1_x = concs
		r1 = -a*x + kd1*a1_x
		dadt = 0.0 +r1 
		da1_xdt = 0.0 -r1 
		dxdt = r1
		return [dadt, dxdt, da1_xdt]
	res = solve_ivp(lambda t, y: ode_multisite_1_to_1(y,t,kd1), interval, [a, x,0.0]).y[2:,-1]
	return (0+1*res[0])/x

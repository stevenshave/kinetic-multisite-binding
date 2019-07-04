import numpy as np
from scipy.integrate import odeint
def multisite_1_to_3(a,x,kd1,kd2,kd3):
	def ode_multisite_1_to_3(concs, t, kd1,kd2,kd3):
		a,x,a1_x,a2_x,a3_x,a1_2_x,a1_3_x,a2_3_x,a1_2_3_x=concs
		r1 = -a*x + kd1*a1_x
		r2 = -a*x + kd2*a2_x
		r3 = -a*x + kd3*a3_x
		r4 = -a1_x*x + kd2*a1_2_x
		r5 = -a2_x*x + kd1*a1_2_x
		r6 = -a1_x*x + kd3*a1_3_x
		r7 = -a3_x*x + kd1*a1_3_x
		r8 = -a2_x*x + kd3*a2_3_x
		r9 = -a3_x*x + kd2*a2_3_x
		r10 = -a1_2_x*x + kd3*a1_2_3_x
		r11 = -a1_3_x*x + kd2*a1_2_3_x
		r12 = -a2_3_x*x + kd1*a1_2_3_x
		dadt = 0.0 +r1 +r2 +r3 
		da1_xdt = 0.0 -r1 +r4 +r6 
		da2_xdt = 0.0 -r2 +r5 +r8 
		da3_xdt = 0.0 -r3 +r7 +r9 
		da1_2_xdt = 0.0 -r4 -r5 +r10 
		da1_3_xdt = 0.0 -r6 -r7 +r11 
		da2_3_xdt = 0.0 -r8 -r9 +r12 
		da1_2_3_xdt = 0.0 -r10 -r11 -r12 
		dxdt = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12
		return [dadt, dxdt, da1_xdt, da2_xdt, da3_xdt, da1_2_xdt, da1_3_xdt, da2_3_xdt, da1_2_3_xdt]
	res = odeint(ode_multisite_1_to_3, [a, x, 0.0,0.0,0.0,0.0,0.0,0.0,0.0], np.linspace(0, 10, 1000), args=(kd1,kd2,kd3))[-1]
	return (1*(sum(res[2:5]))+2*(sum(res[5:8]))+3*res[8])/x

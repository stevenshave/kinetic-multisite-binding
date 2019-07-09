import numpy as np
from scipy.integrate import solve_ivp
def multisite_1_to_4(a, x, kd1, kd2, kd3, kd4, interval=(0,1)):
	def ode_multisite_1_to_4(concs, t, kd1, kd2, kd3, kd4):
		a, x, a1_x, a2_x, a3_x, a4_x, a1_2_x, a1_3_x, a1_4_x, a2_3_x, a2_4_x, a3_4_x, a1_2_3_x, a1_2_4_x, a1_3_4_x, a2_3_4_x, a1_2_3_4_x = concs
		r1 = -a*x + kd1*a1_x
		r2 = -a*x + kd2*a2_x
		r3 = -a*x + kd3*a3_x
		r4 = -a*x + kd4*a4_x
		r5 = -a1_x*x + kd2*a1_2_x
		r6 = -a2_x*x + kd1*a1_2_x
		r7 = -a1_x*x + kd3*a1_3_x
		r8 = -a3_x*x + kd1*a1_3_x
		r9 = -a1_x*x + kd4*a1_4_x
		r10 = -a4_x*x + kd1*a1_4_x
		r11 = -a2_x*x + kd3*a2_3_x
		r12 = -a3_x*x + kd2*a2_3_x
		r13 = -a2_x*x + kd4*a2_4_x
		r14 = -a4_x*x + kd2*a2_4_x
		r15 = -a3_x*x + kd4*a3_4_x
		r16 = -a4_x*x + kd3*a3_4_x
		r17 = -a1_2_x*x + kd3*a1_2_3_x
		r18 = -a1_3_x*x + kd2*a1_2_3_x
		r19 = -a2_3_x*x + kd1*a1_2_3_x
		r20 = -a1_2_x*x + kd4*a1_2_4_x
		r21 = -a1_4_x*x + kd2*a1_2_4_x
		r22 = -a2_4_x*x + kd1*a1_2_4_x
		r23 = -a1_3_x*x + kd4*a1_3_4_x
		r24 = -a1_4_x*x + kd3*a1_3_4_x
		r25 = -a3_4_x*x + kd1*a1_3_4_x
		r26 = -a2_3_x*x + kd4*a2_3_4_x
		r27 = -a2_4_x*x + kd3*a2_3_4_x
		r28 = -a3_4_x*x + kd2*a2_3_4_x
		r29 = -a1_2_3_x*x + kd4*a1_2_3_4_x
		r30 = -a1_2_4_x*x + kd3*a1_2_3_4_x
		r31 = -a1_3_4_x*x + kd2*a1_2_3_4_x
		r32 = -a2_3_4_x*x + kd1*a1_2_3_4_x
		dadt = 0.0 +r1 +r2 +r3 +r4 
		da1_xdt = 0.0 -r1 +r5 +r7 +r9 
		da2_xdt = 0.0 -r2 +r6 +r11 +r13 
		da3_xdt = 0.0 -r3 +r8 +r12 +r15 
		da4_xdt = 0.0 -r4 +r10 +r14 +r16 
		da1_2_xdt = 0.0 -r5 -r6 +r17 +r20 
		da1_3_xdt = 0.0 -r7 -r8 +r18 +r23 
		da1_4_xdt = 0.0 -r9 -r10 +r21 +r24 
		da2_3_xdt = 0.0 -r11 -r12 +r19 +r26 
		da2_4_xdt = 0.0 -r13 -r14 +r22 +r27 
		da3_4_xdt = 0.0 -r15 -r16 +r25 +r28 
		da1_2_3_xdt = 0.0 -r17 -r18 -r19 +r29 
		da1_2_4_xdt = 0.0 -r20 -r21 -r22 +r30 
		da1_3_4_xdt = 0.0 -r23 -r24 -r25 +r31 
		da2_3_4_xdt = 0.0 -r26 -r27 -r28 +r32 
		da1_2_3_4_xdt = 0.0 -r29 -r30 -r31 -r32 
		dxdt = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16+r17+r18+r19+r20+r21+r22+r23+r24+r25+r26+r27+r28+r29+r30+r31+r32
		return [dadt, dxdt, da1_xdt, da2_xdt, da3_xdt, da4_xdt, da1_2_xdt, da1_3_xdt, da1_4_xdt, da2_3_xdt, da2_4_xdt, da3_4_xdt, da1_2_3_xdt, da1_2_4_xdt, da1_3_4_xdt, da2_3_4_xdt, da1_2_3_4_xdt]
	res = solve_ivp(lambda t, y: ode_multisite_1_to_4(y,t,kd1, kd2, kd3, kd4), interval, [a, x,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]).y[2:,-1]
	return (1*(sum(res[0:4]))+2*(sum(res[4:10]))+3*(sum(res[10:14]))+4*res[14])/x

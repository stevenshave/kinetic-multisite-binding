# %%
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
%matplotlib inline
import matplotlib.pyplot as plt


# %%
KON = 1
KD = 5
KOFF = KD/KON


def binding_one_to_one(concs, t, KD):

    KON = 1
    KOFF = KD/KON
    dPdt = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
    dLdt = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
    dPLdt = KON*concs[0]*concs[1]-(KOFF)*concs[2]
    return [dPdt, dLdt, dPLdt]


t = np.linspace(0, 1, 2000)
y0 = [40, 10, 0]
y = odeint(binding_one_to_one, y0, t, args=(5,))
plt.plot(t, y[:, 2]/y0[1])
print(y[:, 2]/y0[1])
# %%


# %%
print("Hello")


def one_to_one_ax(a, x, kdax):
    def ode_one_to_one(concs, t, kdax):
        KON = 1
        KOFF = kdax/KON
        p = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        l = -KON*concs[0]*concs[1]+(KOFF)*concs[2]
        pl = KON*concs[0]*concs[1]-(KOFF)*concs[2]
        return [p, l, pl]

    return odeint(ode_one_to_one, [a, x, 0], np.linspace(0, 1, 100), args=(kdax,))[-1, 2]/x


print(one_to_one_ax(40, 10, 5))


# %%
a = np.linspace(0, 100)
ax = np.empty(len(a))
for i in range(len(a)):
    ax[i] = one_to_one_ax(a[i], 10, 5)
plt.plot(a, ax)
plt.show()
# %%


def one_to_two_ax(a, x, kd1, kd2):
    def ode_one_to_two(concs, t, kda1, kda2):
        a, x, a1x, a2x, a12x = concs
        KON1 = 1.0
        KON2 = 1.0
        KOFF1 = kd1/KON1
        KOFF2 = kd2/KON2

        r1 = -KON1*a*x + KOFF1*a1x
        r2 = -KON2*a*x + KOFF2*a2x
        r3 = -KON2*a1x*x + KOFF2*a12x
        r4 = -KON1*a2x*x + KOFF1*a12x

        dadt = r1+r2
        dxdt = r1+r2+r3+r4
        da1xdt = r3-r1
        da2xdt = r4-r2
        da12xdt = -r3-r4

        return [dadt, dxdt, da1xdt, da2xdt, da12xdt]

    res = odeint(ode_one_to_two, [a, x, 0.0, 0.0, 0.0], np.linspace(
        0, 10, 1000), args=(kd1, kd2))[-1]
    return (res[2]+res[3]+2*res[4])/x


print(one_to_two_ax(40, 10, 5, 5))


# %%
def one_to_three_ax(a, x, kd1, kd2, kd3):
    def ode_one_to_two(concs, t, kda1, kda2, kd3):
        a, x, a1x, a2x, a3x, a12x, a13x, a23x, a123x = concs
        KON1 = 1.0
        KON2 = 1.0
        KON3 = 1.0

        KOFF1 = kd1/KON1
        KOFF2 = kd2/KON2
        KOFF3 = kd3/KON3

        r1 = -a*x + KOFF1*a1x
        r2 = -a*x + KOFF2*a2x
        r3 = -a*x + KOFF3*a3x
        # --
        r4 = -a1x*x + KOFF2*a12x
        r5 = -a1x*x + KOFF3*a13x
        # --
        r6 = -a2x*x+KOFF1*a12x
        r7 = -a2x*x+KOFF3*a13x
        # --
        r8 = -a3x*x+KOFF1*a13x
        r9 = -a3x*x+KOFF2*a23x
        # --
        r10 = -a12x*x+KOFF3*a123x
        r11 = -a13x*x+KOFF2*a123x
        r12 = -a23x*x+KOFF1*a123x

        d_a = r1+r2+r3
        d_x = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12
        d_a1x = -r1+r4+r5
        d_a2x = -r2+r6+r7
        d_a3x = -r3+r8+r9
        d_a12x = -r4-r6+r10
        d_a13x = -r5-r8+r11
        d_a23x = -r7-r9+r12
        d_a123x = -r10-r11-r12
        #        0    1    2     3     4      5      6      7       8
        return [d_a, d_x, d_a1x, d_a2x, d_a3x, d_a12x, d_a13x, d_a23x, d_a123x]

    res = odeint(ode_one_to_two, [a, x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], np.linspace(
        0, 10, 1000), args=(kd1, kd2, kd3))[-1]
    return (res[2] + res[3] + res[4] + 2*(res[5]+res[6]+res[7]) + 3*res[8])/x


print(one_to_three_ax(40, 10, 5, 5,5))


# %%


def one_to_four_ax(a, x, kd1, kd2, kd3, kd4):
    def ode_one_to_four(concs, t, kd1, kd2, kd3, kd4):
        a, x, a1x, a2x, a3x, a4x, a12x, a13x, a14x, a23x, a24x, a34x, a123x, a124x, a134x, a234x, a1234x = concs

        r1 = -a*x + kd1*a1x
        r2 = -a*x + kd2*a2x
        r3 = -a*x + kd3*a3x
        r4 = -a*x + kd4*a4x
        # --
        r5 = -a1x*x + kd2*a12x
        r6 = -a1x*x + kd3*a13x
        r7 = -a1x*x + kd4*a14x
        # --
        r8 = -a2x*x + kd1*a12x
        r9 = -a2x*x + kd3*a23x
        r10 = -a2x*x + kd4*a24x
        # --
        r11 = -a3x*x + kd1*a13x
        r12 = -a3x*x + kd2*a23x
        r13 = -a3x*x + kd4*a34x
        # --
        r14 = -a4x*x + kd1*a14x
        r15 = -a4x*x + kd2*a24x
        r16 = -a4x*x + kd3*a34x
        # --
        r17 = -a12x*x + kd3*a123x
        r18 = -a12x*x + kd4*a124x
        # --
        r19 = -a13x*x + kd2*a123x
        r20 = -a13x*x + kd4*a134x
        # --
        r21 = -a14x*x + kd2*a124x
        r22 = -a14x*x + kd3*a134x
        # --
        r23 = -a23x*x + kd1*a123x
        r24 = -a23x*x + kd4*a234x
        # --
        r25 = -a24x*x + kd1*a124x
        r26 = -a24x*x + kd3*a234x
        # --
        r27 = -a34x*x + kd1*a134x
        r28 = -a34x*x + kd2*a234x
        # --
        r29 = -a123x*x + kd4*a1234x
        r30 = -a124x*x + kd3*a1234x
        r31 = -a134x*x + kd2*a1234x
        r32 = -a234x*x + kd1*a1234x
        
        d_a=r1+r2+r3+r4
        d_x=r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16+r17+r18+r19+r20+r21+r22+r23+r24+r25+r26+r27+r28+r29+r30+r31+r32
        d_a1x=-r1+r5+r6+r7
        d_a2x=-r2+r8+r9+r10
        d_a3x=-r3+r11+r12+r13
        d_a4x=-r4+r14+r15+r16
        d_a12x=-r5-r8+r17+r18
        d_a13x=-r6-r11+r19+r20
        d_a14x=-r7-r14+r21+r22
        d_a23x=-r9-r12+r23+r24
        d_a24x=-r10-r15+r25+r26
        d_a34x=-r13-r16+r27+r28
        d_a123x=-r17-r19-r23+r29
        d_a124x=-r18-r21-r25+r30
        d_a134x=-r20-r22-r27+r31
        d_a234x=-r24-r26-r28+r32
        d_a1234x=-r29-r30-r31-r32
        #       0  1    2   3     4    5    6    7     8      9    10     11   12     13     14      15     16    
        return [d_a, d_x, d_a1x, d_a2x, d_a3x, d_a4x, d_a12x, d_a13x, d_a14x, d_a23x, d_a24x, d_a34x, d_a123x, d_a124x, d_a134x, d_a234x, d_a1234x]

    res = odeint(ode_one_to_four, [a, x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], np.linspace(
        0, 1, 100), args=(kd1, kd2, kd3,kd4))[-1]
    return (sum(res[2:6])+ 2*sum(res[6:12])+ 3*sum(res[12:16]) +4*res[16])/x
print("HERE",one_to_four_ax(40.0, 10.0, 5.0, 5.0,5.0, 5.0))


# %%
a = np.linspace(0, 100, 100)
a1x = np.empty(len(a))
a2x = np.empty(len(a))
a3x = np.empty(len(a))
a4x = np.empty(len(a))
for i in range(len(a)):
    a1x[i] = one_to_one_ax(a[i], 10, 50)
    a2x[i] = one_to_two_ax(a[i], 10, 50, 50)
    a3x[i] = one_to_three_ax(a[i], 10, 50, 50, 50)
    a4x[i] = one_to_four_ax(a[i], 10, 50, 10, 2, 0.1)
fig, ax=plt.subplots(1,1)
ax.plot(a, a4x, 'c',label="4 sites")
ax.plot(a, a3x, 'b',label="3 sites")
ax.plot(a, a2x, 'g',label="2 sites")
ax.plot(a, a1x, 'r',label="1 site")
plt.legend()
plt.xlabel("A")
plt.ylabel("Fraction X bound")


# %%

# %%
from pathlib import Path
from multiprocessing import Pool
import sys
sys.setrecursionlimit(10000)
import numpy as np
import matplotlib.pyplot as plt
import multisite_1_to_1
import multisite_1_to_2
import multisite_1_to_3
import multisite_1_to_4
import multisite_1_to_5
import multisite_1_to_6
import multisite_1_to_7
import multisite_1_to_8
import multisite_1_to_9
import multisite_1_to_10

import pickle
x_conc=40.0
a_conc = np.linspace(0, 100, 200)
kd_for_plot=40
curves = np.zeros((10,len(a_conc)))


class GenParallelPoints:
    def __init__(self, curves, a_conc, xconc, kd_for_plot):
        self.curves=curves
        self.a_conc=a_conc
        self.x_conc=x_conc
        self.kd=kd_for_plot
    def __call__(self, i):
        res=np.zeros(10)
        res[0] = multisite_1_to_1.multisite_1_to_1(self.a_conc[i], x_conc, *([self.kd]*1))
        res[1] = multisite_1_to_2.multisite_1_to_2(self.a_conc[i], x_conc, *([self.kd]*2))
        res[2] = multisite_1_to_3.multisite_1_to_3(self.a_conc[i], x_conc, *([self.kd]*3))
        res[3] = multisite_1_to_4.multisite_1_to_4(self.a_conc[i], x_conc, *([self.kd]*4))
        res[4] = multisite_1_to_5.multisite_1_to_5(self.a_conc[i], x_conc, *([self.kd]*5))
        res[5] = multisite_1_to_6.multisite_1_to_6(self.a_conc[i], x_conc, *([self.kd]*6))
        res[6] = multisite_1_to_7.multisite_1_to_7(self.a_conc[i], x_conc, *([self.kd]*7))
        res[7] = multisite_1_to_8.multisite_1_to_8(self.a_conc[i], x_conc, *([self.kd]*8))
        res[8] = multisite_1_to_9.multisite_1_to_9(self.a_conc[i], x_conc, *([self.kd]*9))
        res[9] = multisite_1_to_10.multisite_1_to_10(self.a_conc[i], x_conc, *([self.kd]*10))
        print(i)
        return res

        
if Path("curves.np").is_file():
    curves=np.load(open("curves.npz", "rb"))
else:
    pool=Pool()
    curves=np.array(pool.map(GenParallelPoints(curves, a_conc, x_conc, kd_for_plot), [i for i in range(len(a_conc))])).T
    np.save(open("curves.npz", "wb"), curves)


#pickle.dump(curves, open("curves.pkl", "wb"))
fig, ax=plt.subplots(1,1)
ax.plot(a_conc, curves[9], label="10 sites")
ax.plot(a_conc, curves[8], label="9 sites")
ax.plot(a_conc, curves[7], label="8 sites")
ax.plot(a_conc, curves[6], label="7 sites")
ax.plot(a_conc, curves[5], label="6 sites")
ax.plot(a_conc, curves[4], label="5 sites")
ax.plot(a_conc, curves[3], label="4 sites")
ax.plot(a_conc, curves[2], label="3 sites")
ax.plot(a_conc, curves[1], label="2 sites")
ax.plot(a_conc, curves[0], label="1 sites")
plt.grid()
plt.legend()
plt.title("Protein-ligand 1:many complexation with multiple ligand binding sites.  40 uM ligand with KDs=40 uM")
plt.xlabel("Protein (uM)")
plt.ylabel("Fraction ligand bound")
ax.set_xlim(a_conc[0], a_conc[-1])
ax.set_ylim(0, 1)

plt.show()
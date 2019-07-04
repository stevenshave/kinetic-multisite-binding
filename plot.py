# %%
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

a = np.linspace(0, 100, 100)
a1x = np.empty(len(a))
a2x = np.empty(len(a))
a3x = np.empty(len(a))
a4x = np.empty(len(a))
a5x = np.empty(len(a))
a6x = np.empty(len(a))
a7x = np.empty(len(a))
a8x = np.empty(len(a))
a9x = np.empty(len(a))
a10x = np.empty(len(a))

kd_for_plot=40.0
for i in range(len(a)):
    a1x[i] = multisite_1_to_1.multisite_1_to_1(a[i], 10, *([kd_for_plot]*1))
    a2x[i] = multisite_1_to_2.multisite_1_to_2(a[i], 10, *([kd_for_plot]*2))
    a3x[i] = multisite_1_to_3.multisite_1_to_3(a[i], 10, *([kd_for_plot]*3))
    a4x[i] = multisite_1_to_4.multisite_1_to_4(a[i], 10, *([kd_for_plot]*4))
    a5x[i] = multisite_1_to_5.multisite_1_to_5(a[i], 10, *([kd_for_plot]*5))
    a6x[i] = multisite_1_to_6.multisite_1_to_6(a[i], 10, *([kd_for_plot]*6))
    a7x[i] = multisite_1_to_7.multisite_1_to_7(a[i], 10, *([kd_for_plot]*7))
    a8x[i] = multisite_1_to_8.multisite_1_to_8(a[i], 10, *([kd_for_plot]*8))
    a9x[i] = multisite_1_to_9.multisite_1_to_9(a[i], 10, *([kd_for_plot]*9))
    a10x[i] = multisite_1_to_10.multisite_1_to_10(a[i], 10, *([kd_for_plot]*10))
    print(i)
pickle.dump(a1x, open("a1x.pkl", "wb"))
pickle.dump(a2x, open("a2x.pkl", "wb"))
pickle.dump(a3x, open("a3x.pkl", "wb"))
pickle.dump(a4x, open("a4x.pkl", "wb"))
pickle.dump(a5x, open("a5x.pkl", "wb"))
pickle.dump(a6x, open("a6x.pkl", "wb"))
pickle.dump(a7x, open("a7x.pkl", "wb"))
pickle.dump(a8x, open("a8x.pkl", "wb"))
pickle.dump(a9x, open("a9x.pkl", "wb"))
pickle.dump(a10x, open("a10x.pkl", "wb"))
fig, ax=plt.subplots(1,1)
ax.plot(a, a10x, label="10 sites")
ax.plot(a, a9x, label="9 sites")
ax.plot(a, a8x, label="8 sites")
ax.plot(a, a7x, label="7 sites")
ax.plot(a, a6x, label="6 sites")
ax.plot(a, a5x, label="5 sites")
ax.plot(a, a4x, label="4 sites")
ax.plot(a, a3x, label="3 sites")
ax.plot(a, a2x, label="2 sites")
ax.plot(a, a1x, label="1 site")
plt.legend()
plt.xlabel("A")
plt.ylabel("Fraction X bound")
plt.show()
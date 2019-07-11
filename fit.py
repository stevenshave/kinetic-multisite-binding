# %%
import lmfit

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

protein_conc=1.0


def _residual(params, protein_conc, x, y, ymin):
    assert len(x)==len(y), "X and Y are different sizes"
    res=0
    
    y_hat=np.zeros(len(y))
    for i in range(len(y_hat)):
        fraction_x_bound=multisite_1_to_3.multisite_1_to_3(protein_conc,x[i], params['kd1'].value,params['kd2'].value,params['kd3'].value)
        amount_x_bound=fraction_x_bound*x[i]
        fraction_total_possible_bound=amount_x_bound/(3*protein_conc)
        y_hat[i]=ymin+((params['ymax'].value-ymin)/protein_conc)*fraction_total_possible_bound
    return np.square(y_hat-y)


x=[0.000010,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2,2.5,5,10,15,20,30]
y=[0.000001,1.01E+08,63134669,45748112,25162485,1.52E+08,2.9E+08,4.3E+08,4.1E+08,7.21E+08,1.04E+09,1.62E+09,1.6E+09,2.1E+09,2.17E+09,2.28E+09]

params=lmfit.Parameters()
params.add("kd1", value=12.0, min=0, max=np.inf)
params.add("kd2", value=20.0, min=0, max=np.inf)
params.add("kd3", value=30.0, min=0, max=np.inf)
params.add("ymax", value=np.max(y), min=0.8*np.max(y), max=2*np.max(y))

lmmini=lmfit.Minimizer(_residual, params, fcn_args=(protein_conc, x,y,np.min(y)))
print(params)

result=lmmini.minimize()
print(result.params)


print("Here")
x_hat=np.linspace(0.00001, np.max(x),num=100)
y_hats=np.zeros(len(x_hat))
for i in range(y_hats.shape[0]):
    
    fraction_x_bound=multisite_1_to_3.multisite_1_to_3(protein_conc,x_hat[i], result.params['kd1'].value,result.params['kd2'].value,result.params['kd3'].value)
    amount_x_bound=fraction_x_bound*x_hat[i]
    fraction_total_possible_bound=amount_x_bound/(3*protein_conc)
    y_hats[i]=np.min(y)+((result.params['ymax'].value-np.min(y))/protein_conc)*fraction_total_possible_bound

print(x_hat)
print(y_hats)

fig, ax=plt.subplots(1,1)
ax.scatter(x,y, 20, marker="+", color='k')
ax.plot(x_hat, y_hats, ':', color='k')
plt.show()


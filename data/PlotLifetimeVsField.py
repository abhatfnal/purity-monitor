import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
import pandas as pd 
from scipy.optimize import curve_fit
import numpy as np
from scipy.interpolate import UnivariateSpline

params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

data = pd.read_csv('/home/aj487/purity/data/20190830.txt', sep=',',dtype=float)
data.columns = ['Voltage', 'Drift Time', 'Collection', 'Lifetime', 'Drift Time Error', 'Collection Error', 'Lifetime Error Low', 'Lifetime Error High']



exo = pd.read_csv('exo200_lifetime_vs_field.txt', sep=',')
exo.columns = ['Electric Field', 'Lifetime', 'Lifetime Low', 'Lifetime High']
ExoLifeTimeError = [exo['Lifetime']-exo['Lifetime Low'], exo['Lifetime High']-exo['Lifetime']]

bakale = pd.read_csv('bakale_lifetime_vs_field.txt', sep=',')
bakale.columns = ['Electric Field', 'Attachment']
BakaleConcentration = 3.6E-9 # mol/l
BakaleLifetime = 1/(bakale['Attachment']*BakaleConcentration) * 1E6

DriftLength = 80 # mm
Resistance = 80 # MOhm
DriftField = data['Voltage'] # total V across DriftLength 

DriftField = (DriftField*Resistance/(200.0+Resistance))/(DriftLength/10.0)


newdata = [x for _,x in sorted(zip(DriftField,data['Lifetime']))]
newfield = sorted(DriftField)
newerror = [x for _,x in sorted(zip(DriftField,data['Lifetime Error Low']))]

nums = [0,1,4,5,8,-1]
newfield = [newfield[i] for i in nums] 
newdata = [newdata[i] for i in nums] 
newerror = np.array([newerror[i] for i in nums])
print(newdata)
spl = UnivariateSpline(newfield, newdata, w=1/newerror, s=5, k=1)
print(DriftField, data['Lifetime'])
X = np.linspace(0,400,1000)

fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.xlim(0,400)
plt.ylim(1E2,5E4)
plt.yscale('log')
plt.xlabel('Drift Field [V/cm]', fontsize=16)
plt.ylabel(r'Lifetime [$\mu$s]', fontsize=16)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.grid(b=True, which='major', color='grey', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
plt.errorbar(DriftField, data['Lifetime'], yerr=[data['Lifetime Error Low'],data['Lifetime Error High']], label='YLXPS', color=colors[0], capsize=3, fmt='o', ms=4)
plt.errorbar(exo['Electric Field'], exo['Lifetime'], yerr=ExoLifeTimeError, label='EXO-200', color=colors[1], capsize=3, fmt='o', ms=4)
plt.plot(bakale['Electric Field'], BakaleLifetime, label='Bakale et. al.', color=colors[2])
plt.plot(X, spl(X), label='Extrapolation', color=colors[0], linestyle='--')

ax.legend(loc='upper right')
fig.tight_layout()
plt.savefig('lifetime_vs_driftfield.pdf')
plt.show()
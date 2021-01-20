import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
import pandas as pd 
from scipy.optimize import curve_fit
import numpy as np

params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

def ExpNeg(x,a,b):
    return 100*np.exp(-(x-a)/b)

data = pd.read_csv('/home/fas/david_moore/aj487/purity_monitor/data/20190830.txt', sep=",")
data.columns = ["Voltage", "Drift Time", "Collection", "Lifetime", 'Drift Time Error', 'Collection Error', 'Lifetime Error Low', 'Lifetime Error High']

DriftLength = 80 # mm
Resistance = 80 # MOhm
DriftTime = data['Drift Time'] # us
DriftField = data['Voltage'] # total V across DriftLength 
DriftField = (DriftField*Resistance/(200.0+Resistance))/(DriftLength/10.0)

fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.xlim(0,400)
plt.ylim(60,120)
plt.xlabel('Drift Field [V/cm]', fontsize=16)
plt.ylabel(r'Charge Collection [\%]', fontsize=16)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.grid(b=True, which='major', color='grey', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
plt.errorbar(DriftField, data['Collection']*100, yerr=data['Collection Error']*100, label='YLXPS', color=colors[0], capsize=2, fmt='o', ms=4)
ax.legend(loc='upper left')
fig.tight_layout()
plt.savefig('collection_vs_driftfield.pdf')


fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.xlim(0,140)
plt.ylim(60,120)
plt.xlabel(r'Drift Time [$\mu$s]', fontsize=16)
plt.ylabel(r'Charge Collection [\%]', fontsize=16)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.grid(b=True, which='major', color='grey', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')


popt, cov = curve_fit(ExpNeg, DriftTime[-4:], data['Collection'][-4:]*100, p0=[100,100])
FitX = np.linspace(72,125,100)
FitY = ExpNeg(FitX, *popt)
plt.errorbar(DriftTime, data['Collection']*100, xerr=[3]*len(DriftTime),yerr=data['Collection Error']*100, label='YLXPS', color=colors[0], capsize=2, fmt='o', ms=4)
plt.plot(FitX, FitY, label=r'$f(t)=100\cdot\exp(-\frac{t-%.1f}{%.1f})$' % (popt[0],popt[1]), color='red')

ax.legend(loc='upper left')
fig.tight_layout()
plt.savefig('collection_vs_drifttime.pdf')

plt.show()
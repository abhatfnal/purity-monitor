import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
import pandas as pd 
from scipy.optimize import curve_fit
import numpy as np

params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

data = pd.read_csv('20190830.txt', sep=",")
data.columns = ["Voltage", "Drift Time", "Collection", "Lifetime", 'Drift Time Error', 'Collection Error', 'Lifetime Error Low', 'Lifetime Error High']

exo = pd.read_csv('exo200_driftvelocity.txt', sep=",")
exo.columns = ["Electric Field", "Drift Velocity"]

DriftLength = 80 # mm
Resistance = 80 # MOhm
DriftTime = data['Drift Time'] # us
DriftField = data['Voltage'] # total V across DriftLength 
DriftField = (DriftField*Resistance/(200.0+Resistance))/(DriftLength/10.0)
DriftVelocity = DriftLength/DriftTime

data['Drift Time Error'] = [3] * len(data['Drift Time Error'])
DriftVelocityError = DriftLength/pow(DriftTime,2)*data['Drift Time Error']

fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.xlim(0,400)
plt.ylim(0,2)
plt.xlabel('Drift Field [V/cm]', fontsize=16)
plt.ylabel(r'Drift Velocity [mm/$\mu$s]', fontsize=16)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.grid(b=True, which='major', color='grey', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
plt.errorbar(DriftField, DriftVelocity, yerr=DriftVelocityError, label='YLXPS', color=colors[0], capsize=2, fmt='o', ms=4)
plt.errorbar(exo['Electric Field'], exo['Drift Velocity'], yerr=0.01*exo['Drift Velocity'], label='EXO-200', color=colors[1], capsize=2, fmt='o', ms=4)
ax.legend(loc='upper left')
fig.tight_layout()
plt.savefig('driftvelocity_vs_driftfield.pdf')
plt.show()
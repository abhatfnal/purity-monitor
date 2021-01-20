# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd 

# params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
# plt.rcParams.update(params) 
# colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

data = pd.read_csv('xe_pressure_vs_density.txt', sep=",")
data.columns = ['Density', 'Pressure']

def func(x,a,b,c,d,e,f):
    # return a*x**3 + b*x**2 + c*x**1 + d 
    return a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x**1 + f

# fig = plt.figure(figsize=(10,6))
# ax = fig.gca()
# plt.xlim(0,2)
# plt.ylim(0,120)
# plt.xlabel(r'Density [g/cm$^3$]', fontsize=16)
# plt.ylabel('Pressure [bar]', fontsize=16)
# ax.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.xaxis.set_major_locator(MultipleLocator(0.2))
# ax.yaxis.set_minor_locator(AutoMinorLocator(5))
# ax.yaxis.set_major_locator(MultipleLocator(10))
# ax.grid(b=True, which='major', color='grey', linestyle='--')
# plt.plot(data['Density'], data['Pressure'], color=colors[0])

popt, pcov = curve_fit(func, data['Density'], data['Pressure'])
print(popt, pcov)
xfit = np.linspace(0.01,0.3,1000)
for x in xfit: 
    print(x, func(x,*popt)*14.5038-14.5038)
quit()
plt.plot(xfit,func(xfit, *popt), ls='--', label='fit', color=colors[1])
ax2 = ax.twinx()
ax2.set_ylabel('Pressure [PSIG]')
ax2.set_ylim(-14.5038,120*14.5038 - 14.5038)
ax2.grid(b=True, which='major', color='grey', linestyle=':')
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_major_locator(MultipleLocator(100))

fig.tight_layout()
plt.savefig('XePressureVsDensity.pdf')
# plt.show()
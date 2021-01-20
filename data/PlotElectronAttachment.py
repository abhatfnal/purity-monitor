import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd 
from scipy.interpolate import interp1d

params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

data = pd.read_csv('electron_attachment.txt', sep=",")
data.columns = ['Electric Field', 'Attachment']

fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.xlim(1E1,1E4)
plt.ylim(1E10,1E12)
plt.xlabel('Electric Field [$\mathrm{V} \mathrm{cm}^{-1}$]', fontsize=14)
plt.ylabel(r'Electron Attachment [$\mathrm{M}^{-1} \mathrm{S}^{-1}$]', fontsize=14)
plt.xscale('log')
plt.yscale('log')
ax.grid(b=True, which='major', color='grey', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')

f2 = interp1d(data['Electric Field'], data['Attachment'], kind='linear')
xnew = np.linspace(np.min(data['Electric Field']), np.max(data['Electric Field']), 100)

for x,y in zip(xnew, f2(xnew)):
    print('%.2f,%.2f' % (x,y))

plt.scatter(data['Electric Field'], data['Attachment'], label='Oxygen', color=colors[0])
plt.plot(xnew, f2(xnew), label='Interpolation', color=colors[1])
ax.legend(loc='upper right')

fig.tight_layout()
plt.savefig('ElectronAttachment.pdf')
plt.show()
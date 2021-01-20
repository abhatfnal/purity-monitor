import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator


plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex':True,'font.size':16,'font.family':'lmodern','text.latex.unicode':True,}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black', 'grey']


fig = plt.figure(0, figsize=(10,6))
ax = fig.gca()
plt.xlabel('Time [s]', fontsize=16)
plt.ylabel('Pressure [mbar]', fontsize=16)
plt.xlim(1E1,1E13)
plt.ylim(1E-10, 1E1)
plt.xscale('log')
plt.yscale('log')
plt.grid()

x1 = np.linspace(1,1E6,1E6)
y1 = 11*np.exp(-0.025*x1)
y2 = 0.1*x1**-1
y3 = 0.001*x1**-0.5
y4 = [1E-9]*len(x1)
plt.plot(x1, y1, color='k')
plt.plot(x1, y2, color='k')
plt.plot(x1, y3, color='k')
plt.plot(x1, y4, color='k')
plt.plot(x1, (y1+y2+y3+y4)*1.1, color='k', linestyle='--')


x2 = np.linspace(1E6,1E13,1E6)
y1 = 11*np.exp(-0.025*x2)
y2 = 0.1*x2**-1
y3 = 0.001*x2**-0.5
y4 = [1E-9]*len(x2)
plt.plot(x2, y1, color='k')
plt.plot(x2, y2, color='k')
plt.plot(x2, y3, color='k')
plt.plot(x2, y4, color='k')
plt.plot(x2, (y1+y2+y3+y4)*1.1, color='k', linestyle='--')

plt.text(0.2, 0.8, 'Volume Gas', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.4, 0.2, 'Desorption', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.63, 0.3, 'Diffusion', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
plt.text(0.25, 0.07, 'Permeation', horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)


# ax.legend(loc='upper right', title='Phase %d - Run %d' % (Phase,jj), prop={'size': 14})
fig.tight_layout()
plt.savefig("pumping.pdf")
plt.close()

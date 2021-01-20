import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex':True,'font.size':16,'font.family':'lmodern','text.latex.unicode':True,}
plt.rcParams.update(params) 

fig = plt.figure(0, figsize=(10,6))
ax = fig.gca()
plt.xlabel('Time [hours]', fontsize=16)
plt.ylabel(r'Outgassing Rate [$\mathrm{mbar} \cdot \mathrm{liter} \cdot \mathrm{s}^{-1} \cdot  \mathrm{cm}^{-2}$]', fontsize=16)
plt.xlim(1E-1,1E2)
plt.ylim(1E-11, 1E-6)
plt.xscale('log')
plt.yscale('log')
plt.grid()

x01 = [0.29, 29]
y01 = [5.5e-8, 1.9e-10]

x02 = [0.29, 15]
y02 = [3.2e-8, 3.9e-10]

x03 = [0.50, 39]
y03 = [1.0e-8, 7.6e-11]

x04 = [0.49, 39]
y04 = [1.0e-8, 6.1e-11]

x05 = [0.5, 39]
y05 = [1.0e-8,  4.1e-11]

x06 = [0.29, 29]
y06 = [ 1.3e-8, 6.4e-11]

x07 = [0.29, 30]
y07 = [3.6e-9, 1.7e-11]

x08 = [0.29, 29]
y08 = [3.1e-9,  1.4e-11]

x09 = [0.2, 20]
y09 = [2.5e-8,  2.5e-10]

x10 = [0.20, 20]
y10 = [1.4e-8, 1.5e-10]

x12 = [0.49, 59]
y12 = [5.4e-8, 1.5e-9]

x13 = [0.49, 59]
y13 = [5.3e-8, 6.0e-10]

x14 = [0.20, 19]
y14 = [4.6e-8, 2.5e-10]

x15 = [0.49, 59]
y15 = [1.1e-8, 2.3e-10]

x16 = [0.49, 59]
y16 = [5.2e-9, 1.4e-10]

plt.plot(x01, y01, color='k')
plt.plot(x02, y02, color='k')
plt.plot(x03, y03, color='k')
plt.plot(x04, y04, color='k')
plt.plot(x05, y05, color='k')
plt.plot(x06, y06, color='k')
plt.plot(x07, y07, color='k')
plt.plot(x08, y08, color='k')
plt.plot(x09, y09, color='k')
plt.plot(x10, y10, color='k')
plt.plot(x12, y12, color='k')
plt.plot(x13, y13, color='k')
plt.plot(x14, y14, color='k')
plt.plot(x15, y15, color='k')

fig.tight_layout()
plt.savefig("desorption.pdf")
plt.close()

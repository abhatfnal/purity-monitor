import numpy as np
import matplotlib
import os, time
# matplotlib.use("PDF")
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
arg = parser.parse_args() 

colors = ['red', 'blue', 'green', 'black', 'orange']

fig = plt.figure(figsize=(12,9))
ax = fig.gca()
ax.grid(b=True, which='major', color='k', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(1))

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel('Mass [u]', fontsize=14)
plt.ylabel('Pressure [Torr]', fontsize=14)
plt.yscale('log')

for jj, x in enumerate(arg.filepath): 
    mass = []
    pressure = []
    file = open(x, "r")
    path, filename = os.path.split(x)
    for ii, line in enumerate(file):
        if(ii<22):
            continue
        else: 
            mass.append(float(line[1:6]))
            pressure.append(float(line[8:-4]))
    plt.plot(mass, pressure, color=colors[jj], label=filename)
    plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)

mass = np.array(mass)
pressure = np.array(pressure)
element = [r'$H_20$', r'$O_2$', r'$N_2$', r'$Ar$', r'$CO_2$']
atomic_num = [18, 32, 28, 40, 44]

for ii, x in enumerate(atomic_num):
    peak = pressure[np.where(mass==x)]
    ax.annotate(element[ii], xy=(x, peak*1.1), xytext=(x-2,peak*3), arrowprops=dict(width=0.5, headwidth=4),fontsize=12)

ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=len(arg.filepath) ,borderaxespad=0, frameon=True, fontsize=12)
plt.xlim(0, 50)
plt.ylim(1e-10, 1e-6)     
plt.savefig("plots/rga.pdf", bbox_inches = 'tight', pad_inches = 0.2)
plt.show()
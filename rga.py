import numpy as np
import matplotlib
import os, time, datetime
# matplotlib.use("PDF")
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator

now = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
arg = parser.parse_args() 

colors = ['red', 'blue', 'green', 'black', 'orange']
element = [r'$\mathrm{H}_2\mathrm{O}$', r'$\mathrm{N}_2$', r'$\mathrm{O}_2$', r'$\mathrm{Ar}$', r'$\mathrm{CO}_2$']
name = ['Water', 'Nitrogen', 'Oxygen', 'Argon', 'Carbondioxide']
atomic_num = [18, 28, 32, 40, 44]

fig = plt.figure(figsize=(12,9))
ax = fig.gca()
ax.grid(b=True, which='major', color='k', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(1))

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel('Mass [u]', fontsize=14)
plt.ylabel('Pressure [mbar]', fontsize=14)
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
            pressure.append(float(line[8:-4])*1.33322)
    plt.plot(mass, pressure, color=colors[jj], label=filename)
    plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)

mass = np.array(mass)
pressure = np.array(pressure)


for ii, x in enumerate(atomic_num):
    peak = pressure[np.where(mass==x)]
    print mass[np.where(mass==x)][0], '\t', peak[0], '\t', name[ii]
    ax.annotate(element[ii], xy=(x, peak*1.1), xytext=(x,peak*3), arrowprops=dict(width=0.5, headwidth=4),fontsize=12)


ax.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
plt.xlim(0, 50)
plt.ylim(1e-10, 1e-6)     
# plt.savefig("plots/rga_"+now+".pdf", bbox_inches = 'tight', pad_inches = 0.2)
plt.show()

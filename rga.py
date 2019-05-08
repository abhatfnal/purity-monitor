import numpy as np
import matplotlib
import os, time, datetime
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator
from decimal import Decimal

now = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath", nargs="*")
arg = parser.parse_args() 

def GetPartialPressure(data, name, atomic_num):
    PumpSpeed = 70*60 # 70 L/min as an average for each gas species 
    x1 = 0.01
    x2 = 14.7e-3 # mbar 
    x3 = 2.3 # L
    x4 = 4 # min
    LeakRate = x1/x2 * x3/x4 # mbar L/min
    print " | Leak Rate:   ", LeakRate, "mbar L/min"
    print " | Printing concentration of gases..."
    print " | ", 'M [u]', '\t', 'P [mbar]', '\t', 'Conc [%]', '\t', 'ppb', '\t\t', 'Name'
    for ii, x in enumerate(data): 
        y = PumpSpeed*x 
        y2 = y/LeakRate/200.0
        y_mass = y2* atomic_num[ii]/136.0 * 1e9
        y2 = '%.3E' % Decimal(y2) # force scientific format
        y_mass = '%.2E' % Decimal(y_mass)
        pressure = '%.3E' % Decimal(x)
        print " | ", atomic_num[ii], '\t\t', pressure, '\t', y2, '\t', y_mass, '\t', name[ii]


colors = ['red', 'blue', 'green', 'black', 'orange']
element = [r'$\mathrm{He}$', r'$\mathrm{CH}_4$', r'$\mathrm{H}_2\mathrm{O}$', r'$\mathrm{N}_2$', r'$\mathrm{O}_2$', r'$\mathrm{Ar}$', r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', r'$\mathrm{CO}_2$']
name = ['Helium', 'Methane', 'Water', 'Nitrogen', 'Oxygen', 'Argon', 'Ethanol','Carbondioxide']
atomic_num = [2, 16, 18, 28, 32, 40, 41, 44]

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

Mass = []
Pressure = []
PPressure = []

for jj, x in enumerate(arg.filepath): 
    mass = []
    pressure = []
    ppressure = []
    file = open(x, "r")
    print " | Analyzing file:   ", x 
    path, filename = os.path.split(x)
    for ii, line in enumerate(file):
        if(ii<22):
            continue
        else: 
            mass.append(float(line[1:6]))
            pressure.append(float(line[8:-4])*1.33322)
    mass = np.array(mass)
    pressure = np.array(pressure)
    for ii, x in enumerate(atomic_num):
        peak = pressure[np.where(mass==x)][0]
        ppressure.append(peak)
    ppressure = np.array(ppressure)

    plt.plot(mass, pressure, color=colors[jj], label=filename, linewidth=1.5)
    plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)

    Mass.append(mass)
    Pressure.append(pressure)
    PPressure.append(ppressure)

    GetPartialPressure(ppressure, name, atomic_num) 

Mass = np.array(Mass)
Pressure = np.array(Pressure)
PPressure = np.array(PPressure)

for ii, x in enumerate(atomic_num):
    peak = np.max(PPressure[:,ii])
    ax.annotate(element[ii], xy=(x, peak*1.2), xytext=(x, peak*3), arrowprops=dict(facecolor='black', width=0.5, headwidth=4),fontsize=14)

ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
plt.xlim(0, 50)
plt.ylim(1e-10, 1e-6)     
# plt.savefig("plots/rga_"+now+".pdf", bbox_inches = 'tight', pad_inches = 0.2)
plt.show()

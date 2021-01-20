import numpy as np
import matplotlib
import os, time, datetime
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator
from decimal import Decimal

now = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
LeakPath = '/home/fas/david_moore/aj487/purity_monitor/plots/leaktest'
SaveName = LeakPath+"/leak_test_"+now+".pdf"

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath", nargs="*")
parser.add_argument("-s", action="store_true", dest="save", default=False)
arg = parser.parse_args() 

if __name__ == '__main__':
    fig = plt.figure(figsize=(12,9))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(10))

    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Time [s]', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    plt.yscale('log')

    for jj, x in enumerate(arg.filepath): 
        time = []
        pressure = []
        file = open(x, "r")
        print(" | Analyzing file:   ", x )
        path, filename = os.path.split(x)
        for ii, line in enumerate(file):
            if(ii<23):
                continue
            else: 
                data = line.split(',')
                # print(ii,data)
                if float(data[1])*1.33322 > 0: 
                    time.append(float(data[0]))
                    pressure.append(float(data[1])*1.33322)
                # time.append(float(line[0:4]))*1.33322)
                # pressure.append(float(line[8:])*1.33322) # convert torr to mbar 
        time = np.array(time)
        pressure = np.array(pressure)

        plt.plot(time[::10], pressure[::10], label=filename, linewidth=1.5)

    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
    plt.xlim(0, np.max(time))  
    plt.ylim(1E-12,1E-8)
    if(arg.save): 
        plt.savefig(SaveName, bbox_inches='tight', pad_inches=0.2)
    plt.show()

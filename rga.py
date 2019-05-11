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
parser.add_argument("-s", action="store_true", dest="save", default=False)
arg = parser.parse_args() 

def GetPartialPressure(data, Dict):
    LeakRate = 0.3911 # mbar L/min
    LXeDensity = 3.1 # g/ml 
    GXeDensity = 0.005761 # g/ml
    XeDensityRatio = LXeDensity/GXeDensity
    CellVolumeRatio = 0.8/2.3 
    print " | Leak Rate:   ", LeakRate, "mbar L/min"
    print " | Printing concentration of gases..."
    print " | ", 'M [u]', '\t', 'P [mbar]', '\t', 'Conc [%]', '\t', 'ppb', '\t\t', 'Name'
    for ii,x in enumerate(list(Dict)):
        Dict[x]['PumpSpeed'] = Dict[x]['PumpSpeed']*60 #convert to min
        y = data[ii] * Dict[x]['PumpSpeed']
        y = y/LeakRate/(XeDensityRatio*CellVolumeRatio)
        y2 = y*Dict[x]['Mass']/136.0 * 1e9
        print " | ", Dict[x]['Mass'], '\t\t', '%.2E'%Decimal(data[ii]), '\t', '%.2E'%Decimal(y), '\t', '%.2E'%Decimal(y2), '\t', x

def GetDictionary():
    #Pump speed is based on https://shop.edwardsvacuum.com/Viewers/Document.ashx?id=2129&lcid=2057
    Dict = {    
        'Helium': {'Mass':4, 'Formula': r'$\mathrm{He}$', 'PumpSpeed':78},
        'Methane': {'Mass':16, 'Formula': r'$\mathrm{CH}_4$', 'PumpSpeed':75.5},
        'Water': {'Mass':18, 'Formula': r'$\mathrm{H}_2\mathrm{O}$', 'PumpSpeed':60},
        'Nitrogen': {'Mass':28, 'Formula': r'$\mathrm{N}_2$', 'PumpSpeed':84},
        'Oxygen': {'Mass':32, 'Formula': r'$\mathrm{O}_2$', 'PumpSpeed':75.5},
        'Argon': {'Mass':40, 'Formula': r'$\mathrm{Ar}$', 'PumpSpeed':80},
        'Ethanol': {'Mass':41, 'Formula': r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', 'PumpSpeed':75.5},
        'Carbondioxide': {'Mass':44, 'Formula': r'$\mathrm{CO}_2$', 'PumpSpeed':75.5}}
    return Dict

if __name__ == '__main__':
    Dict = GetDictionary()

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

    PPressure = []
    for jj, x in enumerate(arg.filepath): 
        mass = []
        pressure = []
        ppressure = []
        
        print " | Analyzing file:   ", x 
        file = open(x, "r")
        path, filename = os.path.split(x)
        for ii, line in enumerate(file):
            if(ii<22): 
                continue
            else: 
                mass.append(float(line[1:6]))
                pressure.append(float(line[8:-4])*1.33322) #torr to mbar
        mass = np.array(mass)
        pressure = np.array(pressure)

        for x in list(Dict):
            peak = pressure[np.where(mass==Dict[x]['Mass'])][0]
            ppressure.append(peak)
        ppressure = np.array(ppressure)
        PPressure.append(ppressure)
        GetPartialPressure(ppressure, Dict) 

        plt.plot(mass, pressure, label=filename, linewidth=1.5)
        plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)

    PPressure = np.array(PPressure)
    for ii, x in enumerate(list(Dict)):
        peak = np.max(PPressure[:,ii])
        ax.annotate(Dict[x]['Formula'], xy=(Dict[x]['Mass'],peak*1.2), xytext=(Dict[x]['Mass'],peak*3), arrowprops=dict(facecolor='black',width=0.5,headwidth=4), fontsize=14)

    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
    plt.xlim(0, 50)
    plt.ylim(1e-10, 1e-6)    
    if(arg.save):  
        plt.savefig("plots/rga_"+now+".pdf", bbox_inches = 'tight', pad_inches = 0.2)
    plt.show()

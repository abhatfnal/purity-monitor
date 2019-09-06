import numpy as np
import os, time, datetime, glob, csv, argparse, h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from decimal import Decimal

import Dictionary as Di
import Plots as Pl

now = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
RGAPath = '/home/fas/david_moore/sjb93'
SaveName = RGAPath+"/rga_"+now+".pdf" 
DataSavePath = '/home/fas/david_moore/aj487/project'
DataSaveName = DataSavePath+"/oxygenpp_"+now+".csv"

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath")
parser.add_argument("-n", type=int, action="store",  dest="num")
arg = parser.parse_args() 

def GetPartialPressure(data, Dict):
    # LeakRate = 0.3911 # mbar L/min
    LeakRate = 115.78 # mbar L/min
    LXeDensity = 3.1 # g/ml 
    GXeDensity = 0.005761 # g/ml
    XeDensityRatio = LXeDensity/GXeDensity
    CellVolumeRatio = 0.8/2.3 
    print(" | Leak Rate:   ", LeakRate, "mbar L/min")
    print(" | Printing concentration of gases...")
    print(" | ", 'M [u]', '\t', 'P [mbar]', '\t', 'Conc [%]', '\t', 'ppb', '\t\t', 'Name')
    for ii,x in enumerate(Dict.keys()):
        Dict[x]['PumpSpeed'] = Dict[x]['PumpSpeed']*60 #convert to min
        y = data[ii]*Dict[x]['PumpSpeed']
        y = y/LeakRate/(XeDensityRatio*CellVolumeRatio)
        y2 = y*Dict[x]['Mass']/136.0*1e9
        print(" | ", Dict[x]['Mass'], '\t', '%.2E'%Decimal(data[ii]), '\t', '%.2E'%Decimal(y), '\t', '%.2E'%Decimal(y2), '\t', x)

def turn_file_names_into_date_times(imaginary_file_names):
    path, imaginary_file_names = os.path.split(imaginary_file_names)
    if imaginary_file_names[-6:-4] == 'PM' or imaginary_file_names[-6:-4] == 'AM':
        return (datetime.datetime.strptime(imaginary_file_names[0:24], '%b_%d_%Y__%I-%M-%S_%p'))
    else: # This is used when the files are given in a 24 hour format
        return (datetime.datetime.strptime(imaginary_file_names[0:21], '%b_%d_%Y__%H-%M-%S'))

def GetData(File): 
    Mass, Pressure, Time = [], [], []

    f = h5py.File(arg.filepath, 'r')
    Keys = list(f.keys())

    Mass = np.array(f.get('Mass'))
    Group = f.get('Pressure')
    GroupKeys = Group.keys()
    print("Reading in %s" % (File))
    print("Number of files %d" % (len(GroupKeys)))
    for ii, key in enumerate(GroupKeys):
        Time.append(key)
        Pressure.append(np.array(Group.get(key))*1.33322)
    return np.array(Mass), np.array(Pressure), np.array(Time)

def GetArrowPosition(Dict, Mass, Pressure): 
    Data = {}
    ArrowPos = {}
    for y in Dict.keys():
        Data[y] = [] 
    for x in Pressure: 
        for y in Dict.keys(): 
            Condition = np.where((Mass > Dict[y]['Mass']-0.5) & (Mass < Dict[y]['Mass']+0.5))
            Values = x[Condition]
            Max = np.max(Values)
            Data[y].append(Max)
    for y in Data.keys(): 
        ArrowPos[y] = np.max(Data[y])
    return Data, ArrowPos


if __name__ == '__main__':
    Dict = Di.GetDictionary()
    LeakRate = Di.GetTurnToLeakRate()
    Mass, Pressure, Time = GetData(arg.filepath)
    Data, ArrowPos = GetArrowPosition(Dict, Mass, Pressure)

    Pl.PlotBestFitOverTime(Time, Data)
    plt.savefig('rga_fit.pdf')
    # plt.show()

    Pl.PlotOverTime(Time, Data)
    plt.savefig('rga_vs_time.pdf')
    # plt.show()

    Pl.PlotRGASpectrum(Time, Mass, Pressure, Dict, ArrowPos)
    plt.savefig('rga_spectrum.pdf')
    # plt.show()

    print("Time elapsed: %.2f sec" % time.process_time())
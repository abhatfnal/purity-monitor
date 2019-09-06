import numpy as np
import matplotlib
import os, time, datetime, glob
import matplotlib.pyplot as plt
import csv
import argparse
from matplotlib.ticker import MultipleLocator
from decimal import Decimal
from matplotlib.dates import DateFormatter
from scipy.optimize import curve_fit
from scipy.stats import chisquare

# Save all measurements in an hdf5 file to speed everything up

now = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
RGAPath = '/home/fas/david_moore/sjb93'
SaveName = RGAPath+"/rga_"+now+".pdf" 
DataSavePath = '/home/fas/david_moore/sjb93/project'
DataSaveName = DataSavePath+"/oxygenpp_"+now+".csv"

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath", nargs="*")
parser.add_argument("-s", action="store_true", dest="save", default=False)
parser.add_argument("-d", type=str, action="store",  dest="dirpath", nargs="*")
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

def GetDictionary():
    #Pump speed is based on https://shop.edwardsvacuum.com/Viewers/Document.ashx?id=2129&lcid=2057
    Dict = {    
        'Helium': {'Mass':4.0, 'Formula': r'$\mathrm{He}$', 'PumpSpeed':78},
        'Carbon': {'Mass':12.0, 'Formula': r'$\mathrm{C}$', 'PumpSpeed':75.5},
        'Methane': {'Mass':16.0, 'Formula': r'$\mathrm{CH}_4$', 'PumpSpeed':75.5},
        'Water': {'Mass':18.0, 'Formula': r'$\mathrm{H}_2\mathrm{O}$', 'PumpSpeed':60},
        'Nitrogen': {'Mass':28.0, 'Formula': r'$\mathrm{N}_2$', 'PumpSpeed':84},
        'Oxygen': {'Mass':32.0, 'Formula': r'$\mathrm{O}_2$', 'PumpSpeed':75.5},
        'Argon': {'Mass':40.0, 'Formula': r'$\mathrm{Ar}$', 'PumpSpeed':80},
        'Ethanol': {'Mass':41.0, 'Formula': r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', 'PumpSpeed':75.5},
        'Carbondioxide': {'Mass':44.0, 'Formula': r'$\mathrm{CO}_2$', 'PumpSpeed':75.5}}
    return Dict

def GetTurnToLeakRate():
    LeakRate = {
        0.40: {'LeakRate':0.09, 'Error': 0.01},
        1.40: {'LeakRate':0.39, 'Error': 0.04},
        1.60: {'LeakRate':2.36, 'Error': 0.26},
        1.65: {'LeakRate':3.91, 'Error': 0.43},
        1.70: {'LeakRate':7.20, 'Error': 0.79},
        1.75: {'LeakRate':11.73, 'Error': 1.29},
        1.80: {'LeakRate':28.16, 'Error': 3.10},
        1.85: {'LeakRate':64.93, 'Error': 7.14},
        1.90: {'LeakRate':115.78, 'Error': 12.74}}
    return LeakRate

def turn_file_names_into_date_times(imaginary_file_names):
    path, imaginary_file_names = os.path.split(imaginary_file_names)
    if imaginary_file_names[-6:-4] == 'PM' or imaginary_file_names[-6:-4] == 'AM':
        return (datetime.datetime.strptime(imaginary_file_names[0:24], '%b_%d_%Y__%I-%M-%S_%p'))
    else: # This is used when the files are given in a 24 hour format
        return (datetime.datetime.strptime(imaginary_file_names[0:21], '%b_%d_%Y__%H-%M-%S'))

def ChooseFilesToAnalyze(arg):
    files = []
    if (arg.dirpath != None):
        for dir in arg.dirpath:
            files.append(glob.glob(dir+"*.txt"))
    if (arg.filepath != None):
        if ('*' in arg.filepath):
            files.append(glob.glob(dir+arg.filepath))
        else:
            files.append(arg.filepath) 	
    files = [val for sublist in files for val in sublist]
    # convert all the file names into datetime strings and then sort them
    return sorted(files, key=turn_file_names_into_date_times)

def PlotOverTime(time, data, files):
    fig = plt.figure(figsize=(12.4,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    # ax.xaxis.set_major_locator(MultipleLocator(2))
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_yscale('log')
    plt.gcf().autofmt_xdate()
    print(str(time[0])[-2:])
    if str(time[0])[-2:] == 'AM' or str(time[0])[-2:] == 'PM':
        formatter = DateFormatter('%b_%d_%Y__%I-%M-%S_%p')
    else: 
        formatter = DateFormatter('%b_%d_%Y__%H-%M-%S')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    if str(time[0])[-2:] == 'AM' or str(time[0])[-2:] == 'PM':
        timet = [datetime.datetime.strptime(x, '%b_%d_%Y__%I-%M-%S_%p') for x in time]
    else: 
        timet = [datetime.datetime.strptime(x, '%b_%d_%Y__%H-%M-%S') for x in time]
    plt.xlim(timet[0], timet[-1])
    label = [r'$\mathrm{H}_2\mathrm{O}$', r'$\mathrm{N}_2$', r'$\mathrm{O}_2$', r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', r'$\mathrm{CO}_2$']
    for ii,x in enumerate(data): 
        plt.plot(timet, x, label=label[ii])
    plt.legend(loc='upper left', prop={'size': 14}, numpoints=1)


def func(x, a, b, c):
    return a * np.exp(-1*b * x) + c

def PlotBestFitOverTime(time,data,files):
    fig = plt.figure(figsize=(12.4,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    if str(time[0])[-2:] == 'AM' or str(time[0])[-2:] == 'PM':
        timet = [datetime.datetime.strptime(x, '%b_%d_%Y__%I-%M-%S_%p') for x in time]
        timet = [datetime.timedelta.total_seconds(x-timet[0]) for x in timet]
        timet = np.array(timet)
    else:
        timet = [datetime.datetime.strptime(x, '%b_%d_%Y__%H-%M-%S') for x in time]
        timet = [datetime.timedelta.total_seconds(x-timet[0]) for x in timet]
        timet = np.array(timet)
    plt.xlim(timet[0], timet[-1])
    label = [r'$\mathrm{H}_2\mathrm{O}$', r'$\mathrm{N}_2$', r'$\mathrm{O}_2$',r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', r'$\mathrm{CO}_2$']
    label2 = ["Fit" r'$\mathrm{H}_2\mathrm{O}$', "Fit" r'$\mathrm{N}_2$', "Fit" r'$\mathrm{O}_2$', "Fit" r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', "Fit" r'$\mathrm{CO}_2$']
    label3 = ["water's fit params: ", "Nitrogen's fit params: ","Oxygen's fit params: ", "Ethanol's fit params: ", "Carbon dioxide fit params: "]
    colors = ['red', 'green', 'blue', 'black', 'purple']
    for ii,x in enumerate(data): 
        plt.scatter(timet, x, label=label[ii], color = colors[ii])
    for ii,x in enumerate(data): 
        param, pcov = curve_fit(func, timet, x, p0=(1E-8,1E-8,1E-8), maxfev=1000000)
        xsquared = chisquare(x,func(timet,*param),axis = None)
        print(str(label3[ii]) + str(param) + " chi squared value = " + str(xsquared) + "\n")
        plt.plot(timet, func(timet, *param), label=label2[ii], color = colors[ii])
    plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)


if __name__ == '__main__':
    files = ChooseFilesToAnalyze(arg)
    Dict = GetDictionary()
    LeakRate = GetTurnToLeakRate()


    fig = plt.figure(figsize=(12.4,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Mass [u]', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    plt.yscale('log')

    PPressure = []
    viton18 = [] 
    viton28 = [] 
    viton32 =[]
    viton41 = []
    viton44 = [] 
    time = [] 
    for jj, x in enumerate(files): 
        mass = []
        pressure = []
        ppressure = []
        
        print(" | Analyzing file:   ", x)
        file = open(x, "r")
        path, filename = os.path.split(x)
        #print filename[-15:-7]
        # print datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        #time.append(filename[-15:-7]) # this line only reads out the hour and minute
        if str(filename[-6:-4]) == 'AM' or str(filename[-6:-4]) == 'PM':
            time.append(filename[0:24])
        else:
            time.append(filename[0:21]) # this line reads out the whole time stamp down to the second
        #print(filename[0:21]) # this line serves as a test
        for ii, line in enumerate(file):
            # This is the line that just skips the first 22 lines of header about the infro of when the data were taken
            if(ii<22): 
                continue
            else: 
                mass.append(float(line[1:6]))
                pressure.append(float(line[8:-3])*1.33322) #torr to mbar, I had to change this line's -4 to a -3 for some reason, the formatting of the files changed?
        mass = np.array(mass)
        pressure = np.array(pressure)

        peak18 = pressure[np.where((mass>18-0.5) & (mass<18+0.5))]
        viton18.append(np.max(peak18))
        peak28 = pressure[np.where((mass>28-0.5) & (mass<28+0.5))]
        viton28.append(np.max(peak28))
        peak32 = pressure[np.where((mass>32-0.5) & (mass<32+0.5))]
        viton32.append(np.max(peak32))
        peak41 = pressure[np.where((mass>41-0.5) & (mass<41+0.5))]
        viton41.append(np.max(peak41))
        peak44 = pressure[np.where((mass>44-0.5) & (mass<44+0.5))]
        viton44.append(np.max(peak44))

        for x in Dict.keys():
            peak = pressure[np.where((mass>Dict[x]['Mass']-0.5) & (mass<Dict[x]['Mass']+0.5))]
            ppressure.append(np.max(peak))
        ppressure = np.array(ppressure)
        PPressure.append(ppressure)
        # GetPartialPressure(ppressure, Dict) 

        plt.plot(mass, pressure, label=path[-8:]+"/"+filename, linewidth=1.5)
        plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)
    """
    file = open('viton.txt', "r")
    mass = []
    pressure = [] 
    for ii, line in enumerate(file):
        print ii, float(line[0:10]), float(line[13:])
        mass.append(float(line[0:10]))
        pressure.append(float(line[13:])*1.33322*100) #torr to mbar
    plt.plot(mass, pressure, label='viton', linewidth=1.5, color="red")
    """
    
    PPressure = np.array(PPressure)
    for ii, x in enumerate(Dict.keys()):
        peak = np.max(PPressure[:,ii])
        ax.annotate(Dict[x]['Formula'], xy=(Dict[x]['Mass'],peak*1.2), xytext=(Dict[x]['Mass'],peak*3), arrowprops=dict(facecolor='black',width=0.8,headwidth=5), fontsize=12, ha='center', bbox=dict(boxstyle="round", fc="white", ec="grey"))
        # ax.text(52,1e-7,'text', bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    # plt.title(datetime.datetime.strptime(path[-8:], '%Y%m%d'))
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=7, borderaxespad=0, fontsize=8)
    plt.xlim(0, 50)
    plt.ylim(1e-11, 1e-5)  
    # plt.xticks(np.arange(0, 50, step=1))
    if(arg.save):  
        plt.savefig(SaveName, bbox_inches='tight', pad_inches=0.2)


    f = open(DataSaveName, 'a+')
    writer = csv.writer(f, lineterminator='\n')
    for x in range(len(time)):
        row = [time[x], viton32[x]]
        writer.writerow(row)
    f.close()



    plt.show()
    PlotOverTime(time, [viton18, viton28, viton32, viton41, viton44], files)
    plt.show()
    PlotBestFitOverTime(time, [viton18, viton28, viton32, viton41, viton44], files)
    plt.show()
import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
from optparse import OptionParser
from itertools import chain
import h5py
import datetime
import dropbox
from matplotlib.dates import DateFormatter
import upload
from matplotlib.ticker import MultipleLocator

startTime = datetime.datetime.now()

usage = 'usage: %prog [options] arg1 arg2'
parser = OptionParser(usage=usage)
parser.add_option('-f', '--file', action='store', type='string', dest='filename', help='Specify path/file to be read in.')
parser.add_option('-l', '--level', action='store', type='string', dest='levelfile', default='', help='Specify level path/file to be read in.')
parser.add_option('-d', '--date', action='store', type='string', dest='date', help='Date of measurement.')
parser.add_option('-t', '--time', action='store', type='int', dest='time', default ='0', help='How much of past data to show in minutes.')
(options, args) = parser.parse_args()

class monitor:
    def __init__(self, filename):
        self.Name = ['Gas System', 'Chamber', 'Stainless-steel Cylinder', 'NOT USED', 'LN Dewar 1', 'LN Dewar 2', 'Xenon Pump', 'Flow Meter', 'Back Pump', 'Cold Head', 'Copper Ring', 'TPC Bottom', 'Copper Jacket',  'Time']
        self.GasSystemP = []
        self.ChamberP = []
        self.SSCylinderP = []
        self.LN1Level = []
        self.LN2Level = []
        self.XenonPumpP = []
        self.FlowMeter = []
        self.BackPumpP = []
        self.ColdHeadT = []
        self.CopperRingT = []
        self.CopperJacketT = []
        self.TPCBottomT = []
        self.Time = []
        self.File = h5py.File(filename, 'r')
        self.LFile = ''
        self.dummy = []
        self.PlotTime = 0

    def GetData(self, data):
        print 'getting data'
        self.GasSystemP = data[:,0]
        self.ChamberP = data[:,1]
        self.SSCylinderP = data[:,2]
        self.dummy = data[:,3]
        self.LN1Level = data[:,4]
        self.LN2Level = data[:,5]
        self.XenonPumpP = data[:,6]
        self.FlowMeter = data[:,7]
        self.BackPumpP = data[:,8]
        self.ColdHeadT = data[:,9]
        self.CopperRingT = data[:,10]
        self.CopperJacketT = data[:,11]
        self.TPCBottomT = data[:,12]
        self.Time = [(x - self.RefTime - 18000) for x in data[:,13]]

def GetDateFromInput(input):
    dd = list(input)
    year = int(dd[0]+dd[1]+dd[2]+dd[3])
    month = int(dd[4]+dd[5])
    day = int(dd[6]+dd[7])
    dt = datetime.datetime(year,month,day,00,00,00)
    tmp = datetime.datetime(1904,1,1,0,0)
    at = int((dt - tmp).total_seconds())
    return at, dt

def PltWfm(Xe,time, data, label, xlabel, ylabel, title, ylim1=1, ylim2=1, yticks=0):
    fig = plt.figure(figsize=(12,9))
    if(ylim1!=1 or ylim2!=1):
        plt.ylim(ylim1, ylim2)
    if(yticks!=0):
        plt.yticks(np.arange(ylim1, ylim2, yticks))
    ax = fig.gca()

    ax.minorticks_on()
    ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=2))
    ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(interval=60))
    ax.yaxis.set_major_locator(MultipleLocator(yticks))
    ax.yaxis.set_minor_locator(MultipleLocator(yticks/2))

    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='k', linestyle=':')

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)

    plt.gcf().autofmt_xdate()
    formatter = DateFormatter('%H:%M')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

    for xx, yy in zip(data,label):
        plt.plot(time, xx, label=yy, linewidth=2)
    plt.legend(loc='upper left')

    if(Xe.PlotTime != 0):
        xlim1 = startTime - datetime.timedelta(seconds=60*Xe.PlotTime)
        xlim2 = startTime + datetime.timedelta(seconds=60*Xe.PlotTime/4)
        ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=10))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(interval=2))
    else:
        xlim1 = Xe.DateTime + datetime.timedelta(seconds=3600*0)
        xlim2 = Xe.DateTime + datetime.timedelta(seconds=3600*24)
    plt.xlim(xlim1, xlim2)

    plt.savefig('MonitorPlots/'+title+'.pdf', bbox_inches = 'tight', pad_inches = 0.2)
    plt.close('all')

def PlotParameters(Xe):
    tt = [(Xe.DateTime + datetime.timedelta(seconds=x)) for x in Xe.Time]

    print 'plotting temperature'
    Temp = [Xe.ColdHeadT, Xe.CopperRingT, Xe.TPCBottomT, Xe.CopperJacketT]
    TempLabels = [Xe.Name[9],Xe.Name[10],Xe.Name[11],Xe.Name[12]]
    PltWfm(Xe, time=tt, data=Temp, label=TempLabels, xlabel='Time [hh:mm]', ylabel='Temperature [C]', title = 'temperature', ylim1=-200, ylim2=30, yticks=20)

    print 'plotting ln level'
    LNLevel = [Xe.LN1Level, Xe.LN2Level]
    LNLabels = [Xe.Name[4],Xe.Name[5]]
    PltWfm(Xe, time=tt, data=LNLevel, label=LNLabels, xlabel='Time [hh:mm]', ylabel='Level [inch]',  title = 'lnlevel', ylim1=0, ylim2=15, yticks=2.)

    print 'plotting flow meter'
    PltWfm(Xe, time=tt, data=[Xe.FlowMeter], label=['Flow Meter'], xlabel='Time [hh:mm]', ylabel='Flow [SLPM]',  title = 'flowmeter', ylim1=0, ylim2=6000, yticks=1000)

    print 'plotting stainless steel bottle pressure'
    Pressure = [Xe.SSCylinderP]
    PressureLabels = [Xe.Name[2]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm]', ylabel='Pressure [PSIG]', title = 'pressure_ssb', ylim1=0, ylim2=1000, yticks=100)

    print 'plotting gas system and chamber pressure'
    Pressure = [Xe.GasSystemP, Xe.ChamberP]
    PressureLabels = [Xe.Name[0], Xe.Name[1]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm]', ylabel='Pressure [PSIG]', title = 'pressure_chamber', ylim1=-15, ylim2=15, yticks=5.)

    print 'plotting pump pressure'
    Pressure = [Xe.XenonPumpP, Xe.BackPumpP]
    PressureLabels = [Xe.Name[6], Xe.Name[8]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm]', ylabel='Pressure [PSIG]', title = 'pressure_pump', ylim1=-15, ylim2=5, yticks=5.)

    if not options.levelfile:
        print 'no capacitance meter data provided'
        pass
    else:
        print 'plotting capacitance meter'
        Xe.LFile = h5py.File(options.levelfile, 'r')
        lkey = list(Xe.LFile.keys())[0]
        lvalues = np.array(Xe.LFile[lkey])
        ltime = [(x - Xe.RefTime - 18000) for x in lvalues[:,1]]
        ldata = lvalues[:,0]
        tt = [(Xe.DateTime + datetime.timedelta(seconds=x)) for x in ltime]
        PltWfm(Xe, tt, [ldata], ['LXe Level'], xlabel='Time [hh:mm]', ylabel='Capacitance [pF]', title = 'lxe_level')

    print 'uploading files to dropbox'
    upload.main()
    print 'finished'

if __name__ == '__main__':

    try:
        Xe = monitor(options.filename)
        key = list(Xe.File.keys())[0]
        values = np.array(Xe.File[key])
        Xe.RefTime, Xe.DateTime = GetDateFromInput(options.date)
        Xe.PlotTime = options.time
        # Xe.PlotTime = 30.0
        Xe.GetData(values)
        PlotParameters(Xe)
    except Exception as e:
        print "....ERROR!"
        print(e)
        pass

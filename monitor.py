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

startTime = datetime.datetime.now()

usage = 'usage: %prog [options] arg1 arg2'
parser = OptionParser(usage=usage)
parser.add_option('-f', '--file', action='store', type='string', dest='filename', help='Specify path/file to be read in.')
parser.add_option('-l', '--level', action='store', type='string', dest='levelfile', default='', help='Specify level path/file to be read in.')
parser.add_option('-d', '--date', action='store', type='string', dest='date', help='Date of measurement.')
(options, args) = parser.parse_args()

class monitor:
    def __init__(self, filename):
        self.Name = ['Gas System', 'Chamber', 'Stainless-steel Cylinder', 'NOT USED', 'LN Dewar 1', 'LN Dewar 2', 'Xenon Pump', 'Flow Meter', 'Back Pump', 'Cold Head', 'Copper Ring', 'Copper Jacket', 'TPC Bottom', 'Time']
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

def PltWfm(Xe,time, data, label, xlabel, ylabel, title):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=270)
    for xx, yy in zip(data,label):
        plt.plot(time, xx, label=yy)
    plt.legend(loc='upper left')
    plt.gcf().autofmt_xdate()
    formatter = DateFormatter('%H:%M:%S')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    # plt.xlim(startTime - datetime.timedelta(seconds=3600*0.1), startTime + datetime.timedelta(seconds=3600*0.1))
    plt.xlim(Xe.DateTime + datetime.timedelta(seconds=3600*10), (Xe.DateTime + datetime.timedelta(seconds=3600*18)))
    # plt.xlim(Xe.DateTime + datetime.timedelta(seconds=3600*12), startTime + datetime.timedelta(seconds=60))
    # plt.show()
    plt.savefig('MonitorPlots/'+title+'.pdf')
    plt.close('all')

def PlotParameters(Xe):
    tt = [(Xe.DateTime + datetime.timedelta(seconds=x)) for x in Xe.Time]

    print 'plotting temperature'
    Temp = [Xe.ColdHeadT, Xe.CopperRingT, Xe.TPCBottomT, Xe.CopperJacketT]
    TempLabels = [Xe.Name[9],Xe.Name[10],Xe.Name[11],Xe.Name[12]]
    PltWfm(Xe, time=tt, data=Temp, label=TempLabels, xlabel='Time [hh:mm:ss]', ylabel='Temperature [C]', title = 'temperature')

    print 'plotting ln level'
    LNLevel = [Xe.LN1Level, Xe.LN2Level]
    LNLabels = [Xe.Name[4],Xe.Name[5]]
    PltWfm(Xe, time=tt, data=LNLevel, label=LNLabels, xlabel='Time [hh:mm:ss]', ylabel='Level [inch]',  title = 'lnlevel')

    print 'plotting flow meter'
    PltWfm(Xe, time=tt, data=[Xe.FlowMeter], label=['Flow Meter'], xlabel='Time [hh:mm:ss]', ylabel='Flow [SLPM]',  title = 'flowmeter')

    print 'plotting stainless steel bottle pressure'
    Pressure = [Xe.SSCylinderP]
    PressureLabels = [Xe.Name[2]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm:ss]', ylabel='Pressure [PSIG]', title = 'pressure_ssb')

    print 'plotting gas system and chamber pressure'
    Pressure = [Xe.GasSystemP, Xe.ChamberP]
    PressureLabels = [Xe.Name[0], Xe.Name[1]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm:ss]', ylabel='Pressure [PSIG]', title = 'pressure_chamber')

    print 'plotting pump pressure'
    Pressure = [Xe.XenonPumpP, Xe.BackPumpP]
    PressureLabels = [Xe.Name[6], Xe.Name[8]]
    PltWfm(Xe, time=tt, data=Pressure, label=PressureLabels, xlabel='Time [hh:mm:ss]', ylabel='Pressure [PSIG]', title = 'pressure_pump')

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
        PltWfm(Xe, tt, [ldata], ['LXe Level'], xlabel='Time [hh:mm:ss]', ylabel='Capacitance [pF]', title = 'lxe_level')

    print 'uploading files to dropbox'
    upload.main()
    print 'finished'

if __name__ == '__main__':

    try:
        Xe = monitor(options.filename)
        key = list(Xe.File.keys())[0]
        values = np.array(Xe.File[key])
        Xe.RefTime, Xe.DateTime = GetDateFromInput(options.date)
        Xe.GetData(values)
        PlotParameters(Xe)
    except:
        print 'something did not work'
        pass

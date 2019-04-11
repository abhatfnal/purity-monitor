import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py
from optparse import OptionParser
from plots import *
from WaveformClass import *
from HelperClasses import *
import TxtExtension as txt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
parser.add_argument("-d", "--dir", action="store_true", dest="dir", default=False, help="Analyze directory and subdirectories.")
parser.add_argument("--txt", action="store_true", dest="txt", default=False, help="")
parser.add_argument('--add', type=str, nargs = '*', help = 'some ids')
arg = parser.parse_args()

def ImportDataFromHDF5(file, channels):
    f = h5py.File(file, 'r')
    Keys = list(f.keys())
    for ch in channels:
        ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
        ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
        Group = f.get(ch.ChName)
        GroupKeys = Group.keys()
        ch.Files += len(GroupKeys)
        print " | Number of files in ch%d...\t" % ch.ID, ch.Files
        for key in GroupKeys:
            ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale)
            ch.TimeStamp.append(Group.get(key).attrs["TimeStamp"])

def DoAnalysis(channels):
    first = True
    for ii, ch in enumerate(channels):
        print " | Processing data in channel %d..." % (ch.ID)
        ch.GetSampling()
        # ch.SetPolarity()
        ch.SubtractBaseline(state=first)
        ch.GetAverageWfm(state=first)
        # ch.GetAllFFT(state=first)
        ch.RemoveNoise(LowCut=0, HighCut=100E3, Order=3, state=first)
        ch.GetAllMaxima(data=ch.Amp, state=first)
        print " | Extremum: ", ch.GetPeak(ch.MeanAmp.tolist())

def SaveNpy(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

def ChangeTime(data):
    return [datetime.datetime.strptime(x, '%H%M%S') for x in data]

if __name__ == '__main__':

    #Initialize class object for each channel.
    ch1 = WFM(arg.filepath, 1)
    ch2 = WFM(arg.filepath, 2)
    ch1.Pol = 1
    ch2.Pol = -1
    channels = [ch1, ch2]

    #Import data from HDF5 file.
    for file in arg.filepath:
        ImportDataFromHDF5(file, channels)

    #Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    DoAnalysis(channels)

    '''Inside here manual analyses can be done and the plotting of data.'''
    print min(ch1.MeanAmp), max(ch1.MeanAmp)
    print min(ch2.MeanAmp), max(ch2.MeanAmp)

    ratio = [-x/y for x,y in zip(ch1.Max, ch2.Max)]
    PltTime(ch1.TimeStamp, [ratio], ['20190410'], 'Time [HH:MM:SS]', 'Ratio')

    # PltWfm(time=ch2.Time, data=[ch1.MeanAmp], legend=['Anode', 'Cathode'])
    # PltWfm(time=ch1.Time, data=list(ch1.Amp), legend=['Anode', 'Cathode'])
    # time2 = np.load('20190410_pur2.npy')[:,0]
    # PltTime(ch1.TimeStamp,[ratio1],['20190410-1'],'Time [HH:MM:SS]','Ratio',create=False)
    # PltTime(time0,[ratio],['20190410-1'],'Time [HH:MM:SS]','Ratio',create=False)
    # PltTime(time1,[ratio1],['20190410-2'],'Time [HH:MM:SS]','Ratio',create=False)
    # PltTime(ch1.TimeStamp,[ratio2],['20190410-2'],'Time [HH:MM:SS]','Ratio',create=False)
    # PltScatter(ch1.TimeStamp, [ch1.Max] , legend=['Anode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    # ch2.FitFullCurve(data=ch2.MeanAmp, start=-300, end=200, repeats=5, state=True)
    # ch2.FitFullCurveDouble(data=ch2.MeanAmp, start=-200, end=400, repeats=5, state=True)

    plt.show()
    print " | Time elapsed: ", time.clock() , "sec"

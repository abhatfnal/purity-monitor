import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py
from optparse import OptionParser
from plots import *
from WaveformClass import *
from HelperClasses import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
parser.add_argument("-d", "--dir", type=str, action="store",  dest="dirpath", nargs="*", help="Analyze directory and subdirectories.")
parser.add_argument("--txt", action="store_true", dest="txt", default=False, help="")
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
            # print type(np.array(Group.get(key)).flatten() * ch.VScale)
            ch.TimeStamp.append(datetime.datetime.strptime(f.attrs['Date']+Group.get(key).attrs["TimeStamp"], '%Y%m%d%H%M%S'))

def DoAnalysis(channels):
    first = True
    for ii, ch in enumerate(channels):
        print " | Processing data in channel %d..." % (ch.ID)
        ch.GetSampling()
        # ch.SetPolarity()
        ch.SubtractBaseline(state=first)
        # ch.GetAverageWfm(state=first)
        # ch.GetAllFFT(state=first)
        # ch.RemoveNoise(LowCut=3, HighCut=100E3, Order=6, state=first)
        ch.GetAllMaxima(data=ch.Amp, state=first)
        # print " | Extremum: ", ch.GetPeak(ch.MeanAmp.tolist())

def SaveNpy(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

def ChangeTime(data):
    return [datetime.datetime.strptime(x, '%Y%m%d%H%M%S') for x in data]

def ChooseFilesToAnalyze(arg):
    files = []
    if (arg.dirpath != None):
        for dir in arg.dirpath:
            files.append(glob.glob(dir+"*.h5"))
    if (arg.filepath != None):
        files.append(arg.filepath)
    files = [val for sublist in files for val in sublist]
    return files

if __name__ == '__main__':

    ###### Initialize class object for each channel.
    ch1 = WFM(1)
    ch2 = WFM(2)
    ch1.Pol = 1
    ch2.Pol = -1
    channels = [ch1, ch2]

    ###### Import data from HDF5 file.
    files = ChooseFilesToAnalyze(arg)
    for file in files:
        ImportDataFromHDF5(file, channels)

    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    DoAnalysis(channels)

    ###### Inside here manual analyses can be done and the plotting of data
    print " | Plotting data..."
    # ratio = - ch1.Max/ch2.Max
    # PltTime(ch1.TimeStamp, [ratio], ['Ratio'], 'Time [HH:MM:SS]', 'Ratio', [0,1])
    # PltTime(ch2.TimeStamp, [ch1.MaxT,ch2.MaxT], ['Anode','Cathode'], 'Time [HH:MM:SS]', 'Peak Time [$\mu$s]', [0,np.max(ch2.MaxT)+30])
    PltTime(ch2.TimeStamp, [ch1.Max,-ch2.Max], ['Anode','Cathode'], 'Time [HH:MM:SS]', 'Amplitude [mV]', [0,np.max(-ch2.Max)+50])
    # print ch2.TimeStamp[0]
    # print ch2.TimeStamp[-1]
    # td = ch2.TimeStamp[-1] - ch2.TimeStamp[0]
    # print type(td)
    #
    # print (td.seconds//60)%60
    # print td.seconds//60.0

    # PltWfm(ch1.Time, [ch1.Amp[0],ch2.Amp[0]], ['Channel 1','Channel 2'])


    # PltWfm(time=ch1.Time, data=list(ch1.AmpClean), legend=['Anode', 'Cathode'])
    # PltWfm(time=ch2.Time, data=list(ch2.AmpClean), legend=['Anode', 'Cathode'])
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

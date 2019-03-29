import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py
from optparse import OptionParser
from plots import *
from WaveformClass import *
from HelperClasses import *
import TxtExtension as txt

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
parser.add_option("-d", "--dir", action="store_true", dest="dir", default=False, help="Analyze directory and subdirectories.")
parser.add_option("-c", "--compare", action="store", type="string", dest="filepath2", default="", help="Specify path/file to be compared to.")
parser.add_option("--txt", action="store_true", dest="txt", default=False, help="")
(options, args) = parser.parse_args()

def ImportDataFromHDF5(file, channels):
    f = h5py.File(file, 'r')
    Keys = list(f.keys())
    for ch in channels:
        ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
        ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
        Group = f.get(ch.ChName)
        GroupKeys = Group.keys()
        ch.Files = len(GroupKeys)
        print " | Number of files:\t", ch.Files
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
        ch.GetAllMaxima(data=ch.AmpClean, state=first)
        print " | Extremum: ", ch.GetPeak(ch.MeanAmp.tolist())

def SaveNpy(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

if __name__ == '__main__':
    #Initialize channel classes for each channel.
    ch1 = WFM(options.filepath, 1)
    ch2 = WFM(options.filepath, 2)
    ch1.Pol = 1
    ch2.Pol = -1
    channels = [ch1, ch2]

    #Import data from HDF5 file.
    ImportDataFromHDF5(options.filepath, channels)

    #Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    DoAnalysis(channels)

    '''Inside here manual analyses can be done and the plotting of data.'''

    # ratio = [-x/y for x,y in zip(ch1.Max, ch2.Max)]

    # new = detrend(ch2.Amp[0], type='constant', bp=10)
    # PltWfm(time=ch2.Time, data=[ch2.Amp[0]], legend=['Anode', 'Cathode'])
    # PltWfm(time=ch2.Time, data=[new], legend=['Anode', 'Cathode'])
    # PltTime(ch2.TimeStamp,[ratio],'Anode/Cathode','Minutes','Ratio')


    # print min(ch1.MeanAmp), max(ch1.MeanAmp)
    # print min(ch2.MeanAmp), max(ch2.MeanAmp)


    # PltScatter(ch1.TimeStamp, [ch1.Max] , legend=['Anode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    # PltScatter(ch2.TimeStamp, [ch2.Max] , legend=['Cathode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    # PltScatter(ch1.TimeStamp, [ratio] , legend=['Anode/Cathode'], xlabel='Time %H%M%S', ylabel='Ratio')

    # PltTime(ch2.TimeStamp, [ch2.Max], legend='Cathode', xlabel='Minutes', ylabel='Amplitude [mV]')
    # PltScatterD(seg, ch2.Max, ch3.Max, label='Anode vs Time', xlabel='count', ylabel='Amplitude [mV]')
    # PltWfm(time=ch2.Time, data=[ch1.MeanAmp, ch2.MeanAmp], legend=['Anode', 'Cathode'])

    # ch2.FitFullCurve(data=ch2.MeanAmp, start=-300, end=200, repeats=5, state=True)
    # ch2.FitFullCurveDouble(data=ch2.MeanAmp, start=-200, end=400, repeats=5, state=True)

    plt.show()
    print " | Time elapsed: ", time.clock() , "sec"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import time, datetime, sys, os, glob, struct, h5py
from optparse import OptionParser
from itertools import chain
from scipy import signal,integrate
from plots import *
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from WaveformClass import WFM
from HelperClasses import *

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
parser.add_option("-d", "--dir", action="store_true", dest="dir", default=False, help="Analyze directory and subdirectories.")
parser.add_option("-c", "--compare", action="store", type="string", dest="filepath2", default="", help="Specify path/file to be compared to.")
parser.add_option("--txt", action="store_true", dest="txt", default=False, help="")
(options, args) = parser.parse_args()

def ReadFilesInDirectories(filepath):
    directories = sorted(os.listdir(filepath))
    files = []
    NrFiles = 0
    for dir in directories:
        NewDir = filepath+"/"+dir+"/"
        print " | Path:\t", NewDir
        NrFiles2, Files = ReadFilesInDirectory(NewDir)
        files.append(Files)
        NrFiles += NrFiles2
    files = [val for sublist in files for val in sublist]
    return NrFiles, files

def ReadFilesInDirectory(filepath):
    files = glob.glob(filepath+"*.txt")
    files2 = sorted(files)
    NrFiles = len(files)
    print " | Number of files:\t", NrFiles
    return NrFiles, files2

def GetFileCreationTime(file):
    timestamp = os.path.getctime(file)
    modificationTime = time.strftime('%H:%M:%S', time.localtime(timestamp))
    return modificationTime

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

def ReadData(channels, files):
    print " | Reading in data files..."
    for i,line in enumerate(files):
        ProgressBar(i+1, len(files))
        with open(line) as wfm:
            for j,data in enumerate(wfm):
                columns = data.split()
                for k, ch in enumerate(channels):
                    if(j==0):
                        ch.CreateTime.append(GetFileCreationTime(line))
                    if(i==0):
                        ch.Time.append(float(columns[0])*ch.TScale)
                    ch.Amp[i].append(float(columns[ch.ID])*ch.VScale)

def GetNpyData(file):
    data = np.load(file)
    return data

def SubtractTemplate(time, data, file, start, end, state=False):
    if(state): print " | Subtracting template..."
    backg = np.load(file)
    background = backg[:,1]
    DataMinusFit = [(x-y) for x,y in zip(data[start:end],background[start:end])]
    new = data[:start]
    new.extend(DataMinusFit)
    new.extend(data[end:])
    return new

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

def SaveFileWtihNumPyArray(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

if __name__ == '__main__':
    #Initialize channel classes for each channel. Currently one channel (ch1) is the trigger and the other one (ch2) is the signal of both anode and cathode.
    ch1 = WFM(options.filepath, 1)
    ch2 = WFM(options.filepath, 2)
    ch1.Pol = 1
    ch2.Pol = -1

    if(options.txt):
        #Read in files in a given directory. Filename is the file that contains the names of all the data files. num_lines is the number of such data files.
        if(options.dir):
            num_files, files = ReadFilesInDirectories(options.filepath)
        else:
            num_files, files = ReadFilesInDirectory(options.filepath)
        #Get all waveforms in data files and save them in lists. The content of each waveform is saved in chX.Amp[i], with X the channel number and i the waveform number. The time is saved only once in chX.Time
        ReadData([ch1, ch2], files)
        if(options.filepath2 != ""):
            num_files2, files2 = ReadFilesInDirectory(options.filepath2)
            ch3 = WFM(num_files2, options.filepath2, 2)
            ReadData([ch3], files2)
    else:
        #New Method of getting data from HDF5 files.
        ImportDataFromHDF5(options.filepath, [ch1,ch2])

    # print len(ch2.Time)
    # print len(ch1.Amp[0])
    # print len(ch2.Amp[0])

    PltWfm(time=ch2.Time, data=[ch1.Amp[0], ch2.Amp[0]], legend=['Anode', 'Cathode'])

    #Inside here the analysis is carried out, including baseline subtraction, averaging of waveforms, getting the fourier spectrum, applying a frequency bandpass filter and finding the maxima of each waveform.
    if(options.filepath2 != ""):
        DoAnalysis([ch2, ch3])
    else:
        DoAnalysis([ch1, ch2])

    PltWfm(time=ch2.Time, data=[ch1.Amp[0], ch2.Amp[0]], legend=['Anode', 'Cathode'])
    PltWfm(time=ch2.Time, data=[ch1.AmpClean[0], ch2.AmpClean[0]], legend=['Anode', 'Cathode'])


    print min(ch1.MeanAmp), max(ch1.MeanAmp)
    print min(ch2.MeanAmp), max(ch2.MeanAmp)
    ratio = [-x/y for x,y in zip(ch1.Max, ch2.Max)]
    # PltScatter(ch1.TimeStamp, [ch1.Max, ch2.Max] , legend=['Anode', 'Cathode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    # PltScatter(ch1.TimeStamp, [ch1.Max] , legend=['Anode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    # PltScatter(ch2.TimeStamp, [ch2.Max] , legend=['Cathode'], xlabel='Time %H%M%S', ylabel='Amplitude [mV]')
    PltScatter(ch1.TimeStamp, [ratio] , legend=['Anode/Cathode'], xlabel='Time %H%M%S', ylabel='Ratio')
    PltTime(ch1.TimeStamp, [ratio], legend='Anode', xlabel='Minutes', ylabel='Ratio')
    # PltTime(ch2.TimeStamp, [ch2.Max], legend='Cathode', xlabel='Minutes', ylabel='Amplitude [mV]')
    # PltScatterD(seg, ch2.Max, ch3.Max, label='Anode vs Time', xlabel='count', ylabel='Amplitude [mV]')
    # PltWfm(time=ch2.Time, data=[ch1.MeanAmp, ch2.MeanAmp], legend=['Anode', 'Cathode'])

    # ch2.FitFullCurve(data=ch2.MeanAmp, start=-300, end=200, repeats=5, state=True)
    # ch2.FitFullCurveDouble(data=ch2.MeanAmp, start=-200, end=400, repeats=5, state=True)

    plt.show()
    print " | Time elapsed: ", time.clock() , "sec"

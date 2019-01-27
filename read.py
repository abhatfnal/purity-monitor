import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib.pyplot as plt
from optparse import OptionParser
from itertools import chain
from scipy import signal,integrate
from plots import PPltWfm, PltWfm, PltScatter
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from wv_class import WFM
import h5py

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
parser.add_option("-p", "--plot", action="store_true", dest="plot", default=False, help="Show fit plots.")
parser.add_option("-s", "--show", action="store_true", dest="show", default=False, help="Show waveforms.")
parser.add_option("-c", "--compare", action="store", type="string", dest="filepath2", default="", help="Specify path/file to be compared to.")
(options, args) = parser.parse_args()

def ReadFilesInDirectory(filepath):
    filename = "list.dat"
    os.system("ls "+filepath+"*.txt > "+filename)
    num_files = sum(1 for line in open(filename))
    print " | Number of files:          ", num_files
    return num_files, filename

def ReadData(channels, filename):
    file = open(filename, "r")
    print " | Reading in data files..."
    for i,line in enumerate(file):
        # ProgressBar(i+1, ch1.Files)
        wfm = open(line[:-1])
        for j,data in enumerate(wfm):
            columns = data.split()
            for k, ch in enumerate(channels):
                if(i==0):
                    ch.Time.append(float(columns[0])*ch.TScale)
                ch.Amp[i].append(float(columns[ch.ID])*ch.VScale)
        wfm.close()
    file.close()

    for ch in channels:
        ch.GetSampling()

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
        print " | Processing data in channel %d..." % (ii+1)
        ch.SubtractBaseline(state=first)
        ch.GetAverageWfm(state=first)
        # ch.GetAllFFT(state=first)
        # ch.RemoveNoise(LowCut=0, HighCut=100E3, Order=9, state=first)
        ch.GetAllMaxima(data=ch.Amp, state=first)
        ch.Plot = options.plot

def FindTime(data, bin):
    return np.abs((np.asarray(data)-bin)).argmin()

def SaveFileWtihNumPyArray(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

if __name__ == '__main__':
    #Read in files in a given directory. Filename is the file that contains the names of all the data files. num_lines is the number of such data files.
    num_files, filename = ReadFilesInDirectory(options.filepath)

    #Initialize channel classes for each channel. Currently one channel (ch1) is the trigger and the other one (ch2) is the signal of both anode and cathode.
    ch1 = WFM(num_files, 1, "Trigger", 1)
    ch2 = WFM(num_files, 1, "Signal", 2)
    ReadData([ch1,ch2], filename)
    if(options.filepath2 != ""):
        num_files2, filename2 = ReadFilesInDirectory(options.filepath2)
        ch3 = WFM(num_files2, 1, "Anode", 2)
        ReadData([ch3], filename2)
    #Get all waveforms in data files and save them in lists. The content of each waveform is saved in chX.Amp[i], with X the channel number and i the waveform number. The time is saved only once in chX.Time


    #Inside here the analysis is carried out, including baseline subtraction, averaging of waveforms, getting the fourier spectrum, applying a frequency bandpass filter and finding the maxima of each waveform.
    if(options.filepath2 != ""):
        DoAnalysis([ch2, ch3])
    else:
        DoAnalysis([ch2])

    print ch2.Max
    print ch2.MaxT

    PPltWfm(time=ch2.Time, data=ch2.MeanAmp, data2=ch3.MeanAmp , label='purified', label2='normal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')
    # PltScatter(xvalue=ch2.MaxT, yvalue=ch2.Max, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')
    # new = SubtractTemplate(ch2.Time, ch2.MeanAmp, 'noise.npy', FindTime(ch2.Time, -500), FindTime(ch2.Time, 500))
    # PltWfm(time=ch2.Time, data=new, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')

    # PPltWfm(time=ch2.Time, data=new, data2=filt2, label='Cathode', label2='Anode', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')

    # filt = ch2.RemoveNoiseSingle(ch2.MeanAmp, 100, 5E5, 3)
    # filt2 = ch3.RemoveNoiseSingle(ch3.MeanAmp, 100, 5E5, 3)
    #
    # print " | Extrema: ", np.min(ch2.MeanAmp), np.max(ch2.MeanAmp)
    # ch2.FitFullCurveDouble(data=ch2.MeanAmp, start=-450, end=400, repeats=5, state=True)
    # SaveFileWtihNumPyArray(ch2.Time, ch2.MeanAmp, 'noise')

    # _,_,_,_,_ = ch2.FitFullCurve(data=ch2.MeanAmp, start=-200, end=800, repeats=5, state=True)
    # _,_,_,_,_ = ch3.FitFullCurve(data=ch3.MeanAmp, start=-300, end=500, repeats=5, state=True)

    # anode = GetNpyData('anode.npy')
    # anode = anode[:,1].tolist()
    # cathode = GetNpyData('cathode.npy')
    # cathode = cathode[:,1].tolist()

    # PPltWfm(time=ch2.Time, data=ch2.MeanAmp, data2=filt , label='raw', label2='filt', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', save=False)
    # _,_,_,_,_ = ch2.FitFullCurve(data=filt, start=-500, end=500, repeats=5, state=True)


    # ch2.FitFullCurve(data=filt, start=-42, end=-10, repeats=10, state=True)
    # expo = ch2.FitExponential(ch2.MeanAmp, 500, 2500, 3, state=True)
    # test = ch2.SubtractFunction(ch2.MeanAmp, expo, FindTime(ch2.Time, 12), FindTime(ch2.Time, 1000))


    print " | Time elapsed: ", time.clock() , "sec"

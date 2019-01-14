import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib.pyplot as plt
from optparse import OptionParser
from itertools import chain
from scipy import signal,integrate
from plots import PPltWfm, PltWfm
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from wv_class import WFM
import h5py

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
parser.add_option("-p", "--plot", action="store_true", dest="plot", help="Show fit plots.")
parser.add_option("-s", "--show", action="store_true", dest="show", help="Show waveforms.")
(options, args) = parser.parse_args()

def ReadFilesInDirectory():
    filename = "list.dat"
    os.system("ls "+options.filepath+"*.txt > "+filename)
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
                ch.Amp[i].append(float(columns[k+1])*ch.VScale)
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
    print backg
    print background
    DataMinusFit = [(x-y) for x,y in zip(data[start:end],background[start:end])]
    new = data[:start]
    new.extend(DataMinusFit)
    new.extend(data[end:])
    return new

def DoAnalysis(channels):
    first = True
    for ch in channels:
        # PltWfm(time=ch.Time, data=ch.Amp[0][:], label="Signal", xlabel="Time", ylabel="Amp")
        ch.SubtractBaseline(state=first)
        ch.GetAverageWfm(state=first)
        # ch.GetAllFFT(state=first)
        # ch.RemoveNoise(LowCut=0, HighCut=100E3, Order=9, state=first)
        # ch.GetAllMaxima(data=ch.AmpClean, state=first)
        ch.Plot = options.plot
        first = False

def FindTime(data, bin):
    return np.abs((np.asarray(data)-bin)).argmin()

def SaveFileWtihNumPyArray(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

if __name__ == '__main__':
    #Read in files in a given directory. Filename is the file that contains the names of all the data files. num_lines is the number of such data files.
    num_files, filename = ReadFilesInDirectory()

    #Initialize channel classes for each channel. Currently one channel (ch1) is the trigger and the other one (ch2) is the signal of both anode and cathode.
    ch1 = WFM(num_files, -1, "Cathode")
    ch2 = WFM(num_files, 1, "Anode")

    #Get all waveforms in data files and save them in lists. The content of each waveform is saved in chX.Amp[i], with X the channel number and i the waveform number. The time is saved only once in chX.Time
    ReadData([ch1, ch2], filename)

    #Inside here the analysis is carried out, including baseline subtraction, averaging of waveforms, getting the fourier spectrum, applying a frequency bandpass filter and finding the maxima of each waveform.
    DoAnalysis([ch2])



    # PltWfm(time=ch2.Time, data=ch2.MeanAmp, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')

    # new = SubtractTemplate(ch2.Time, ch2.MeanAmp, 'background_model_test.npy', FindTime(ch2.Time, 12), FindTime(ch2.Time, 100))
    # PltWfm(time=ch2.Time, data=new, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]' )

    # SaveFileWtihNumPyArray(ch2.Time, ch2.MeanAmp, 'background_model_test')
    # _,_,_,_,_ = ch2.FitFullCurve(data=ch2.MeanAmp, start=-500, end=1000, repeats=7, state=True)
    # print np.min(ch2.MeanAmp)
    _,_,_,_,_ = ch2.FitFullCurve(data=ch2.MeanAmp, start=-500, end=800, repeats=5, state=True)
    # ch2.FitFullCurveDouble(data=ch2.MeanAmp, start=-450, end=450, repeats=5, state=True)


    quit()

    #Apply a bandpass filter to the first waveform in ch2 with lowcut, highcut and order and then plot it.
    filt = ch2.RemoveNoiseSingle(ch2.MeanAmp, 50000, 5E6, 3)
    PltWfm(time=ch2.Time, data=filt, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', save=True)
    quit()
    # ch2.FitFullCurve(data=filt, start=-42, end=-10, repeats=10, state=True)
    expo = ch2.FitExponential(ch2.MeanAmp, 500, 2500, 3, state=True)

    test = ch2.SubtractFunction(ch2.MeanAmp, expo, FindTime(ch2.Time, 12), FindTime(ch2.Time, 1000))

    PltWfm(time=ch2.Time, data=test, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]')
    # PltWfm(time=ch2.Time, data=filt, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', xlim=1, xlim2=1, ylim=1, ylim2=1)


    print " | Time elapsed: ", time.clock() , "sec"

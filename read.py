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

def ReadData(ch1, ch2, filename):
    file = open(filename, "r")
    print " | Reading in data files..."
    for i,line in enumerate(file):
        # if(i==5):break
        # ProgressBar(i+1, ch1.Files)
        # print i, line[:-1]
        wfm = open(line[:-1])
        for j,data in enumerate(wfm):
            columns = data.split()
            if(i==0):
                ch1.Time.append(float(columns[0])*ch1.TScale)
                ch2.Time.append(float(columns[0])*ch2.TScale)
            ch1.Amp[i].append(float(columns[1])*ch1.VScale)
            ch2.Amp[i].append(float(columns[2])*ch2.VScale)
        wfm.close()
    file.close()
    for ch in (ch1,ch2):
        ch.GetSampling()

def DoAnalysis(ch1, ch2):
    first = False
    for ch in (ch1,ch2):
        # PltWfm(time=ch.Time, data=ch.Amp[0][:], label="Signal", xlabel="Time", ylabel="Amp")
        ch.SubtractBaseline(state=first)
        ch.GetAverageWfm(state=first)
        ch.GetAllFFT(state=first)
        ch.RemoveNoise(LowCut=0.01, HighCut=100E3, Order=9, state=first)
        ch.GetAllMaxima(data=ch.AmpClean, state=first)
        ch.Plot = options.plot
        first = False

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
    ReadData(ch1, ch2, filename)

    #Inside here the analysis is carried out, including baseline subtraction, averaging of waveforms, getting the fourier spectrum, applying a frequency bandpass filter and finding the maxima of each waveform.
    DoAnalysis(ch1, ch2)

    #Apply a bandpass filter to the first waveform in ch2 with lowcut, highcut and order and then plot it.
    filt = ch2.RemoveNoiseSingle(ch2.Amp[0], 0.01, 100E3,12)
    # PltWfm(time=ch2.Time, data=filt, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', xlim=1,xlim2=1,ylim=-10,ylim2=10)

    # filt = ch2.butter_bandpass_filter(ch2.MeanAmp, 500, 100E3, ch2.Sampling, 12).tolist()
    # PltWfm(time=ch2.Time, data=filt, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', xlim=1,xlim2=1,ylim=-10,ylim2=10)
    # filt = ch2.butter_bandpass_filter(ch2.MeanAmp, 1000, 100E3, ch2.Sampling, 12).tolist()
    # PltWfm(time=ch2.Time, data=filt, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]', xlim=1,xlim2=1,ylim=-10,ylim2=10)
    quit()

    SaveFileWtihNumPyArray(self.Time, self.Amp[0])

    new = ch1.RemoveNoiseSingle(ch1.MeanAmp, LowCut=0, HighCut=300E3, Order=3)
    new2 = ch2.RemoveNoiseSingle(ch2.MeanAmp, LowCut=0, HighCut=300E3, Order=3)

    _,_ = ch1.GetPeak(new, ch1.Pol)
    _,_ = ch2.GetPeak(new2, ch2.Pol)

    new3 = ch1.ShapeGaussian(np.asarray(new),sigma=15)      #gaussian convoluted signal
    new4 = ch2.ShapeGaussian(np.asarray(new2),sigma=15)     #gaussian convoluted signal

    start, amp = ch1.FitExponential(new, ch1.PeakTime, ch1.Time[-2000], 2, state=True)
    start2, amp2 = ch2.FitExponential(new2, ch2.PeakTime, ch2.Time[-2000], 7, state=False)

    _,_ = ch1.GetPeak(new3, ch1.Pol)
    _,_ = ch2.GetPeak(new4, ch2.Pol)

    print " | Drift Time Difference:   ", (ch2.PeakTime-ch1.PeakTime) - (start2-start)
    print " | Amplitude Ratio Difference:  ", abs(amp/amp2)-abs(min(new3)/max(new4))
    print " | Time elapsed: ", time.clock() , "sec"
    print " | Field Ratio:  ", abs(ch2.Peak/ch1.Peak)
    print ch1.PeakTime, "\t", ch2.PeakTime, "\t", start, "\t", start2, "\t", ch1.Peak, "\t", ch2.Peak, "\t", amp, "\t", amp2

    if(options.show):
        PPltWfm(ch1.Time,ch1.MeanAmp,ch2.MeanAmp,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)
        PPltWfm(ch1.Time,new,new2,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)
        PPltWfm(ch1.Time,new3,new4,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)

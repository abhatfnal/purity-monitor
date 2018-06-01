import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib.pyplot as plt
from optparse import OptionParser
from itertools import chain
from scipy import signal,integrate
from plots import PPltWfm
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from wv_class import WFM
from features import ProgressBar

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
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
        ProgressBar(i+1, ch1.Files)
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
    for data in (ch1,ch2):
        data.GetSampling()

def DoAnalysis(ch1, ch2):
    first = True
    for data in (ch1,ch2):
        data.SubtractBaseline(state=first)
        data.GetAverageWfm(state=first)
        data.GetAllFFT(state=first)
        data.RemoveNoise(LowCut=0, HighCut=400E3, Order=12, state=first)
        data.GetAllMaxima(data.AmpClean, state=first)
        first = False


if __name__ == '__main__':
    num_files, filename = ReadFilesInDirectory()
    ch1 = WFM(num_files, -1, "Cathode")
    ch2 = WFM(num_files, 1, "Anode")

    ReadData(ch1, ch2, filename)
    DoAnalysis(ch1, ch2)

    new = ch1.RemoveNoiseSingle(ch1.MeanAmp, LowCut=0, HighCut=300E3, Order=3)
    new2 = ch2.RemoveNoiseSingle(ch2.MeanAmp, LowCut=0, HighCut=300E3, Order=3)

    new3 = ch1.ShapeGaussian(np.asarray(new),sigma=15)      #gaussian convoluted signal
    new4 = ch2.ShapeGaussian(np.asarray(new2),sigma=15)     #gaussian convoluted signal

    print " | Fitting exponential to decaying edge..."
    start, amp = ch1.FitExponential(new, ch1.PeakTime, ch1.Time[-1], 2)
    start2, amp2 = ch2.FitExponential(new2, ch2.PeakTime, ch2.Time[-1], 5)

    _,_ = ch1.GetPeak(new3, ch1.Pol)
    _,_ = ch2.GetPeak(new4, ch2.Pol)


    print " | Drift Time Difference:   ", (ch2.PeakTime- ch1.PeakTime) - (start2-start)
    print " | Amplitude Ratio Difference:  ", abs(amp/amp2)-abs(min(new3)/max(new4))
    print ch1.PeakTime, "\t", ch2.PeakTime, "\t", start, "\t", start2, "\t", ch1.Peak, "\t", ch2.Peak, "\t", amp, "\t", amp2
    
    # PPltWfm(ch1.Time,ch1.MeanAmp,ch2.MeanAmp,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)
    # PPltWfm(ch1.Time,new,new2,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)
    # PPltWfm(ch1.Time,new3,new4,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',scale=1.2)

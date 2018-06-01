import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib.pyplot as plt
from optparse import OptionParser
from scipy.signal import blackman
from itertools import chain
from scipy import signal,integrate
from plots import plot_double_waveform, plot_single_waveform, plot_single_fft
from scipy.fftpack import fft
from scipy.optimize import curve_fit
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

    ch1.GetSampling()
    ch2.GetSampling()

if __name__ == '__main__':
    num_files, filename = ReadFilesInDirectory()
    ch1 = WFM(num_files, -1, "Cathode")
    ch2 = WFM(num_files, 1, "Anode")

    ReadData(ch1, ch2, filename)

    print " | Subtracting baseline..."
    ch1.SubtractBaseline()
    ch2.SubtractBaseline()

    print " | Getting average waveform..."
    ch1.GetAverageWfm()
    ch2.GetAverageWfm()

    print " | Getting Fourier spectra..."
    ch1.GetAllFFT()
    ch2.GetAllFFT()

    print " | Removing noise..."
    ch1.RemoveNoise(LowCut=0, HighCut=400E3, Order=12)
    ch2.RemoveNoise(LowCut=0, HighCut=400E3, Order=12)

    print " | Getting extrema of individual files..."
    ch1.GetAllMaxima(ch1.AmpClean)
    ch2.GetAllMaxima(ch2.AmpClean)

    new = ch1.RemoveNoiseSingle(ch1.MeanAmp, LowCut=0, HighCut=300E3, Order=3)
    new2 = ch2.RemoveNoiseSingle(ch2.MeanAmp, LowCut=0, HighCut=300E3, Order=3)

    new3 = ch1.ShapeGaussian(np.asarray(new),sigma=15)      #gaussian convoluted signal
    new4 = ch2.ShapeGaussian(np.asarray(new2),sigma=15)     #gaussian convoluted signal

    print " | Fitting exponential to decaying edge..."
    start, amp = ch1.FitExponential(new, ch1.PeakTime, ch1.Time[-1], 5)
    start2, amp2 = ch2.FitExponential(new2, ch2.PeakTime, ch2.Time[-1], 5)

    _,_ = ch1.GetPeak(new3, ch1.Pol)
    _,_ = ch2.GetPeak(new4, ch2.Pol)

    print ch1.PeakTime, "\t", ch2.PeakTime, "\t", start, "\t", start2, "\t", ch1.Peak, "\t", ch2.Peak, "\t", amp, "\t", amp2
    print " | Drift Time Difference:   ", (ch2.PeakTime- ch1.PeakTime) - (start2-start)
    # print ch1.Time[new3.index(min(new3))], "\t", ch2.Time[new4.index(max(new4))], "\t", start, "\t", start2, "\t", min(new3), "\t", max(new4), "\t", amp, "\t", amp2
    # print " | Drift Time Difference:   ", (ch2.Time[new4.index(max(new4))] - ch1.Time[new3.index(min(new3))]) - (start2-start)
    print " | Amplitude Ratio Difference:  ", abs(amp/amp2)-abs(min(new3)/max(new4))

    plot_double_waveform(ch1.Time,ch1.MeanAmp,ch2.MeanAmp,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',ch1.Time[0],ch1.Time[-1],min(min(ch1.MeanAmp),min(ch2.MeanAmp))*1.2,max(max(ch1.MeanAmp),max(ch2.MeanAmp))*1.2)
    plot_double_waveform(ch1.Time,new,new2,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',ch1.Time[0],ch1.Time[-1],min(min(new),min(new2))*1.2,max(max(new),max(new2))*1.2)
    plot_double_waveform(ch1.Time,new3,new4,'Cathode','Anode','Time [$\mu$s]','Amplitude [mV]',ch1.Time[0],ch1.Time[-1],min(min(new3),min(new4))*1.2,max(max(new3),max(new4))*1.2)


    # for i in range(ch1.Files):
    #     plot_double_waveform(ch1.Time, ch1.AmpClean[i], ch2.AmpClean[i], ch1.label, ch2.label,'Time [$\mu$z]', 'Amplitude [mV]', min(ch1.Time), max(ch1.Time), -10, 20)

    # plot_single_fft(ch1.TimeFFT, ch1.AmpFFT[0], ch2.AmpFFT[0], ch1.label, ch2.label,'Frequency [Hz]', 'Amplitude', min(ch1.TimeFFT), max(ch1.TimeFFT), min(ch1.AmpFFT[0])*2, max(ch1.AmpFFT[0])*2)

    # plot_double_waveform(ch1.Time, ch1.Amp[0], ch2.Amp[0], 'Cathode', 'Anode','Time [$\mu$s]', 'Amplitude [mV]',ch1.Time[0],ch1.Time[-1],-150,50)

    # plot_double_waveform(ch1.Time, new, new2, 'Cathode', 'Anode','Time [$\mu$s]', 'Amplitude [mV]',ch1.Time[0],ch1.Time[-1],-150,50)

    # plot_double_waveform(ch1.Time, ch1.MeanAmp, ch2.MeanAmp, 'Cathode', 'Anode','Time [$\mu$s]', 'Amplitude [mV]',ch1.Time[0],ch1.Time[-1], min(ch1.MeanAmp), max(ch2.MeanAmp))

    # plot_double_waveform(np.linspace(0.0, ch1.Files, float(ch1.Files), endpoint=True), ch1.Max, ch2.Max, ch1.label, ch2.label,'File [#]', 'Amplitude [mV]', 0, ch1.Files, min(ch1.Max)*2, max(ch2.Max)*2)

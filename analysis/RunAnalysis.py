import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py, argparse

from WaveformClass import *
from HelperClasses import *

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
parser.add_argument("-d", "--dir", type=str, action="store",  dest="dirpath", nargs="*", help="Analyze directory and subdirectories.")
parser.add_argument("--txt", action="store_true", dest="txt", default=False, help="")
parser.add_argument("-s", action="store_true", dest="show", default=False, help="")
arg = parser.parse_args()

def ImportDataFromHDF5(file, channels):
    f = h5py.File(file, 'r')
    print " | Filename...", file 
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
            ch.TimeStamp.append(datetime.datetime.strptime(f.attrs['Date']+Group.get(key).attrs["TimeStamp"], '%Y%m%d%H%M%S'))
        
def DoAnalysis(channels):
    first = True
    for ii, ch in enumerate(channels):
        print " | Processing data in channel %d..." % (ch.ID)
        ch.Amp = np.array(ch.Amp)
        ch.GetSampling()
        ch.SubtractBaseline(state=first)
        # ch.GetAllFFT(state=first)
        ch.RemoveNoise(LowCut=0, HighCut=100E3, Order=1, state=first)
        # ch.SubtractLinearBackground(ch.Time, ch.AmpClean, state=first)
        ch.GetAllMaxima(data=ch.AmpClean, state=first)
        # ch.GetIntegral(Data=ch.AmpClean)
        ch.GetAverageWfm(Data=ch.AmpClean, state=first)
    print " | Time elapsed: ", time.clock() , "sec"

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
        if ('*' in arg.filepath):
            files.append(glob.glob(dir+arg.filepath))
        else:
            files.append(arg.filepath)
    files = [val for sublist in files for val in sublist]
    return files

def StandardPlots(ch1, ch2):
    print " | Plotting data..."
    ratio = -ch1.Max/ch2.Max
    DriftTime = ch1.MaxT - ch2.MaxT
    print " | Drift Time...", np.mean(DriftTime), np.std(DriftTime)
    print " | Charge collection...", np.mean(ratio), np.std(ratio), np.std(ratio)/np.sqrt(len(ratio))
    print " | Lifetime...", -np.mean(DriftTime)/np.log(np.mean(ratio)), 'ms'
    print " | Time elapsed: ", time.clock() , "sec"

    DiffMinute = int((np.max(ch1.TimeStamp) - np.min(ch1.TimeStamp)).seconds/60.0 + 0.5)
    XTicks = [int(DiffMinute/12.0/3.0*60.0), int(DiffMinute/12.0 + 0.5)]

    PltTime(Time=ch1.TimeStamp,Data=[ch1.Max,-1*ch2.Max,ratio*100],Legend=['Anode','Cathode','Charge Collection [\%]'],Label='Amplitude [mV]',XTicks=XTicks,YTicks=[5,20],Save='amp_ratio')
    PltTime(Time=ch1.TimeStamp,Data=[ch1.MaxT,ch2.MaxT,DriftTime],Legend=['Anode','Cathode','Drift Time [$\mu$s]'],Label='Peak Time [$\mu$s]',XTicks=XTicks,YTicks=[10,50],Save='drift_time')
    PltTime(Time=ch1.TimeStamp,Data=[ch1.BaseStd,ch2.BaseStd],Legend=['Anode','Cathode'],Label='Baseline Noise [mV]',XTicks=XTicks,YTicks=[1,1],YRange=[0,4],Save='baseline')
    PltTime(Time=ch1.TimeStamp,Data=[ratio],Legend=[''],Label='Charge Collection',YRange=[0,2],XTicks=XTicks,YTicks=[0.1,0.5],Save='ratio')
    PltWfm(Time=ch1.Time,Data=[ch1.MeanAmp,-1*ch2.MeanAmp],Legend=['Anode','Cathode'],XTicks=[50,100],YTicks=[5,10],Save='avg_waveform')
    # PltWfm(Time=ch1.Time,Data=list(ch1.Amp),Legend=[''],XTicks=[50,100],YTicks=[5,10],Color='k',Save='ch1_all_waveforms')
    # PltWfm(Time=ch2.Time,Data=list(ch2.Amp),Legend=[''],XTicks=[50,100],YTicks=[5,10],Color='k',Save='ch2_all_waveforms')

    if arg.show: 
        plt.show()
    else: 
        plt.close('all')


if __name__ == '__main__':

    ###### Initialize class object for each channel.
    ch1 = WFM(1)
    ch2 = WFM(2)
    ch1.Pol = +1
    ch2.Pol = -1
    channels = [ch1, ch2]

    ###### Import data from HDF5 file.
    files = ChooseFilesToAnalyze(arg)
    for file in files:
        ImportDataFromHDF5(file, channels)

    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    DoAnalysis(channels)

    ###### Plotting and saving standard plots for the analysis 
    StandardPlots(ch1,ch2)
    
    print " | Time elapsed: ", time.clock() , "sec"
   

    ###### Additional manual analysis and plotting 


    # for x in ch2.Amp: 

    #     rel = np.min(ch2.MeanAmp)/np.min(x)
    #     corr = signal.correlate(x, ch2.MeanAmp/rel, mode='same') / np.sum(x**2)
    #     fig, (ax_orig, ax_filt, ax_corr) = plt.subplots(3, 1, sharex=True)
    #     ax_orig.plot(x)
    #     ax_orig.set_title('Original signal')
    #     ax_filt.plot(ch2.MeanAmp)
    #     ax_filt.set_title('Template signal')
    #     # ax_corr.plot(corr)
    #     ax_corr.plot(x-ch2.MeanAmp)
    #     ax_corr.set_title('Auto-correlated')
    #     ax_orig.margins(0, 0.1)
    #     fig.tight_layout()
    #     fig.show()
    #     raw_input()

    # quit() 



    # FitFullCurve(ch2, data=ch2.MeanAmp, start=-300, end=500, repeats=5, state=True)

    # PltWfm(Time=ch1.Time,Data=[ch2.Amp[10],ch2.AmpClean[10],ch2.AmpSubtract[10]],Legend=['Raw','Filtered','Bkgd Subtracted'],Label=['Time [HH:MM:SS]','Amplitude [mV]'],YRange=[np.min(ch2.MeanAmp)*1.5,np.max(ch1.MeanAmp)*2.5],XTicks=[20,100],YTicks=[5,20])

    # PltWfm(Time=ch2.Time,Data=list(ch2.Amp),Legend=['Waveforms'],Label=['Time [HH:MM:SS]','Amplitude [mV]'],XTicks=[20,100],YTicks=[5,20],Color='k')
    # PltWfm(Time=ch1.Time,Data=[ch2.MeanAmp],Legend=['Cathode'],Label=['Time [HH:MM:SS]','Amplitude [mV]'],XTicks=[20,100],YTicks=[5,20])
    # plt.show() 
    # XTicks=[300,10]
    # Legend = ['Anode','Cathode','Ratio x100']
    # bins = np.linspace(0.6,1.4,40)
    # ratio = -ch1.Max/ch2.Max
    # plt.hist(ratio, bins, alpha=0.8, label='Raw', color='red')
    # plt.axvline(x=np.mean(ratio), linewidth=2, color='red')
    # test1 = OptimalFilter(ch1.Time, ch1.AmpClean, ch1.Pol)
    # test2 = OptimalFilter(ch2.Time, ch2.AmpClean, ch2.Pol)
    # ratiot = test1/test2
    # plt.hist(ratiot, bins, alpha=0.8, label='Raw', color='blue')
    # plt.axvline(x=np.mean(ratiot), linewidth=2, color='blue')
    # plt.show() 
    # PltTime(Time=ch1.TimeStamp,Data=[test1,test2,ratiot*100],Legend=Legend, Label=['Time [HH:MM:SS]','Amplitude [mV]'],XTicks=XTicks,YTicks=[5,20])
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.Max,-1*ch2.Max,ratio*100],Legend=Legend, Label=['Time [HH:MM:SS]','Amplitude [mV]'],XTicks=XTicks,YTicks=[5,20])
    # plt.show()



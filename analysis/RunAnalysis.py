import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py, argparse

from WaveformClass import *
from HelperClasses import *
from MatchedFilter import MatchedFilter

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
parser.add_argument("-d", "--dir", type=str, action="store",  dest="dirpath", nargs="*", help="Analyze directory and subdirectories.")
parser.add_argument("--txt", action="store_true", dest="txt", default=False, help="")
parser.add_argument("-s", action="store_true", dest="show", default=False, help="")
arg = parser.parse_args()

def ImportDataFromHDF5(file, channels):
    f = h5py.File(file, 'r')
    print(" | Filename...", file)
    Keys = list(f.keys())
    for ch in channels:
        ch.Amp = []
        ch.TimeStamp = []
        ch.Voltages.append(GetVoltageFromName(file))
        ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
        ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
        Group = f.get(ch.ChName)
        GroupKeys = Group.keys()
        ch.Files = len(GroupKeys)
        # print(" | Number of files in ch%d...\t" % ch.ID, ch.Files)
        for key in GroupKeys:
            ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale)
            ch.TimeStamp.append(datetime.datetime.strptime((f.attrs['Date']+Group.get(key).attrs["TimeStamp"]).decode('utf-8'), '%Y%m%d%H%M%S'))
        
def DoAnalysis(channels):
    Print = False
    for ii, ch in enumerate(channels):
        # print(" | Processing data in channel %d..." % (ch.ID))
        ch.Amp = np.array(ch.Amp)
        ch.GetSampling()
        ch.SubtractBaseline(state=Print)
        # ch.GetAllFFT(state=Print)
        ch.RemoveNoise(LowCut=1E2, HighCut=50E3, Order=3, state=Print)
        # ch.SubtractLinearBackground(ch.Time, ch.AmpClean, state=Print)
        ch.GetAllMaxima(data=ch.AmpClean, state=Print)
        # ch.GetIntegral(Data=ch.AmpClean)
        ch.GetAverageWfm(Data=ch.AmpClean, state=Print)
        ch.FindMaxGradient(Data=ch.AmpClean ,state=Print)
    print(" | Time elapsed: ", time.process_time() , "sec")

def GetVoltageFromName(file):
    StartChar = file.find('-')
    EndChar = file.find('V')
    FirstVoltage = int(file[StartChar+1:EndChar])
    NewFile = file[EndChar+1:]
    StartChar = NewFile.find('-')
    EndChar = NewFile.find('V')
    SecondVoltage = int(NewFile[StartChar+1:EndChar])
    TotalVoltage = FirstVoltage + SecondVoltage 
    return TotalVoltage

def SaveNpy(time, data, filename='test'):
    data = np.column_stack((np.asarray(time), np.asarray(data)))
    np.save(filename, data)

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

def PrintResults(ch1, ch2): 
    ChargeCollection = - ch1.Max / ch2.Max
    DriftTime = ch1.GradTime - ch2.GradTime
    
    ChargeCollectionMean = np.mean(ChargeCollection)
    DriftTimeMean = np.mean(DriftTime)
    LifeTime = - np.mean(DriftTime) / np.log(ChargeCollectionMean)

    ChargeCollectionError = np.std(ChargeCollection)/np.sqrt(len(ChargeCollection) )
    DriftTimeError = np.std(DriftTime)
    LifeTimeErr = np.sqrt( pow(DriftTimeError/np.log(ChargeCollectionMean), 2) + pow(DriftTimeMean * ChargeCollectionError / ( pow(np.log(ChargeCollectionMean),2) * ChargeCollectionMean), 2) )   

    print("%.2f,%.3f,%.2f,%.1f,%.3f,%.1f" % (DriftTimeMean, ChargeCollectionMean, LifeTime, DriftTimeError, ChargeCollectionError, LifeTimeErr ) ) 

def StandardPlots(ch1, ch2):
    # print(" | Plotting data...")

    ch1.TimeStamp = np.array(ch1.TimeStamp)
    good =  np.where(ch1.BaseStd<15)

    ch1.Max = ch1.Max[good]
    ch2.Max = ch2.Max[good]

    ch1.MaxT = ch1.MaxT[good]
    ch2.MaxT = ch2.MaxT[good]

    ch1.GradTime = ch1.GradTime[good]
    ch2.GradTime = ch2.GradTime[good]

    ch1.BaseStd = ch1.BaseStd[good]
    ch2.BaseStd = ch2.BaseStd[good]

    ch1.TimeStamp = ch1.TimeStamp[good]

    ratio = -ch1.Max/ch2.Max
    PrintResults(ch1, ch2)
    # DriftTime = ch1.MaxT - ch2.MaxT + 25
    DriftTime = ch1.GradTime-ch2.GradTime
    print(" | Drift Time...", np.mean(DriftTime), np.std(DriftTime))
    print(" | Charge collection...", np.mean(ratio), np.std(ratio), np.std(ratio)/np.sqrt(len(ratio)))
    print(" | Lifetime...", -np.mean(DriftTime)/np.log(np.mean(ratio)), 'us')

    # print("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" % (np.mean(DriftTime), np.mean(ratio), -np.mean(DriftTime)/np.log(np.mean(ratio)), np.std(DriftTime)), np.std(ratio)/np.sqrt(len(ratio)))) )
    # print(" | Drift Time...", np.median(DriftTime), np.std(DriftTime))
    # print(" | Charge collection...", np.median(ratio), np.std(ratio), np.std(ratio)/np.sqrt(len(ratio)))
    # print(" | Lifetime...", -np.median(DriftTime)/np.log(np.median(ratio)), 'us')
    # print(" | Time elapsed: ", time.process_time() , "sec")

    DiffMinute = int((np.max(ch1.TimeStamp) - np.min(ch1.TimeStamp)).seconds/60.0 + 0.5)
    XTicks = int((DiffMinute/12.0 + 0.5))+1

    SavePath = '/home/fas/david_moore/aj487/purity_monitor/plots/analysis/'
    Date = datetime.datetime.now().strftime("%Y%m%d")
    Save = False 
    if Save:
        if not os.path.exists(SavePath+Date):
            os.makedirs(SavePath+Date)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.Max,-1*ch2.Max,ratio*100],Legend=['Anode','Cathode','Charge Collection [\%]'],Label='Amplitude [mV]',XTicks=XTicks,YTicks=50,SaveName='amp_ratio',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.MaxT,ch2.MaxT,DriftTime],Legend=['Anode','Cathode','Drift Time [$\mu$s]'],Label='Peak Time [$\mu$s]',XTicks=XTicks,YTicks=50,YRange=[0,350],SaveName='drift_time',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.GradTime,ch2.GradTime,DriftTime],Legend=['AnodeG','CathodeG','Drift Time G [$\mu$s]'],Label='Peak Time [$\mu$s]',XTicks=XTicks,YTicks=50,YRange=[0,350],SaveName='drift_time',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.BaseStd,ch2.BaseStd],Legend=['Anode','Cathode'],Label='Baseline Noise [mV]',XTicks=XTicks,YTicks=2,YRange=[0,10],SaveName='baseline',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ratio],Legend=[''],Label='Charge Collection',YRange=[0,2],XTicks=XTicks,YTicks=0.5,SaveName='ratio',Save=Save)
    # PltWfm(Time=ch1.Time,Data=[ch1.MeanAmp,-1*ch2.MeanAmp],Legend=['Anode','Cathode'],XTicks=100,YTicks=50,SaveName='avg_waveform')
    # PltWfm(Time=ch1.Time,Data=list(ch1.AmpClean),Legend=['Anode','Cathode'],XTicks=100,YTicks=50,SaveName='avg_waveform',Color='k')
    # PltAllWfm(Time=ch2.Time,Data=list(ch2.Amp)[:10],Legend=[''],XTicks=100,YTicks=50,Save='all_waveform',Color='k')
    
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
        print(file)
        GetVoltageFromName(file)
        ImportDataFromHDF5(file, channels)
        DoAnalysis(channels)
        PrintResults(ch1, ch2)
        # StandardPlots(ch1,ch2)
    quit()

    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    DoAnalysis(channels)

    ###### Plotting and saving standard plots for the analysis 
    PrintResults(ch1, ch2)
    # StandardPlots(ch1,ch2)
    
    print(" | Time elapsed: ", time.process_time() , "sec")
    quit() 

    test = []
    test2 = []
    space = 50
    x = np.linspace(0,len(ch1.Amp[0][::space])*space,len(ch1.Amp[0][::space]))
    for amp in ch1.Amp: 
        prime = np.gradient(amp[::space])
        peak = np.where(prime==np.max(prime))[0][0]
        newx = int(x[peak])
        test.append(ch1.Time[newx])
    for amp in ch2.Amp: 
        prime = np.gradient(amp[::space])
        peak = np.where(prime==np.min(prime))[0][0]
        newx = int(x[peak])
        test2.append(ch2.Time[newx])
    PltTime(Time=ch1.TimeStamp,Data=[ch1.MaxT,ch2.MaxT,test,test2],Legend=['Anode','Cathode','Anode Grad','Cathode Grad'],Label='Peak Time [$\mu$s]',XTicks=1,YTicks=50,YRange=[0,100],SaveName='drift_time',Save=False)
    ch1.MaxT = np.array(ch1.MaxT)
    ch2.MaxT = np.array(ch2.MaxT)
    print(test)
    print(test2)
    quit()

    test = np.array(test)
    test2 = np.array(test2)
    print(np.mean(ch1.MaxT-ch2.MaxT), np.mean(test-test2))

    PltTime(Time=ch1.TimeStamp,Data=[ch1.MaxT-ch2.MaxT, test-test2],Legend=['Drift Time','Drift Time Grad'],Label='Peak Time [$\mu$s]',XTicks=1,YTicks=50,YRange=[0,100],SaveName='drift_time',Save=False)
    plt.show()


    quit()

    ###### Additional manual analysis and plotting 
    # dx = ch1.Sampling 
    # dy = np.diff(ch1.Amp[0])
    space = 50
    prime = np.gradient(ch1.Amp[0][::space], 1)
    prime = prime * np.max(ch1.Amp[0])/np.max(prime)
    plt.plot(np.linspace(0,len(prime)*space,len(prime)), prime, color='red')
    plt.plot(ch1.Amp[0])
    plt.show()
    plt.close()
    prime = np.gradient(-ch2.Amp[0][::space], 1)
    prime = prime * np.max(-ch2.Amp[0])/np.max(prime)
    plt.plot(np.linspace(0,len(prime)*space,len(prime)), prime, color='red')
    plt.plot(-ch2.Amp[0])
    plt.show()
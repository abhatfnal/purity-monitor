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
        ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
        ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
        Group = f.get(ch.ChName)
        GroupKeys = Group.keys()
        ch.Files.append(len(GroupKeys))
        print(" | Number of files in ch%d...\t %d/%d" % (ch.ID, ch.Files[-1], np.sum(ch.Files)))
        for key in GroupKeys:
            ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale * ch.Pol)
            ch.TimeStamp.append(datetime.datetime.strptime((f.attrs['Date']+Group.get(key).attrs["TimeStamp"]).decode('utf-8'), '%Y%m%d%H%M%S'))
        
def DoAnalysis(channels):
    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
    Print = False 
    for ii, ch in enumerate(channels):
        print(" | Processing data in channel %d..." % (ch.ID))
        ch.Amp = np.array(ch.Amp)
        ch.TimeStamp = np.array(ch.TimeStamp)
        ch.GetSampling()
        ch.SubtractBaseline(state=Print)
        ch.RemoveNoise(LowCut=1E1, HighCut=5E4, Order=3, state=Print)
        ch.GetAllMaxima(data=ch.AmpClean, state=Print)
        ch.FindMaxGradient(Data=ch.AmpClean ,state=Print)
        # ch.ApplyCut(Cut=np.where(channels[0].BaseStd<10), state=Print)
        ch.GetAverageWfm(Data=ch.AmpClean, state=Print)
        ch.GetBaselineNoise(Data=ch.AmpClean, state=Print)
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
    FileNames = []
    Method = ''
    if (arg.dirpath != None):
        for dir in arg.dirpath:
            FileNames.append(glob.glob(dir+"*.h5"))
        Method = 'Single'
    if (arg.filepath != None):
        if ('*' in arg.filepath):
            FileNames.append(glob.glob(dir+arg.filepath))
        else:
            FileNames.append(arg.filepath)
        Method = 'All'
    FileNames = [val for sublist in FileNames for val in sublist]
    return FileNames, Method

def PrintResults(ch1, ch2): 
    ChargeCollection = ch1.Max / ch2.Max
    DriftTime = ch1.GradTime - ch2.GradTime
    
    ChargeCollectionMean = np.mean(ChargeCollection)
    DriftTimeMean = np.mean(DriftTime)
    LifeTime = - np.mean(DriftTime) / np.log(ChargeCollectionMean)

    ChargeCollectionError = np.std(ChargeCollection)/np.sqrt(len(ChargeCollection) )
    DriftTimeError = np.std(DriftTime)
    LifeTimeErr = np.sqrt(pow(DriftTimeError/np.log(ChargeCollectionMean),2)+pow(DriftTimeMean*ChargeCollectionError/(pow(np.log(ChargeCollectionMean),2)*ChargeCollectionMean),2))   
    LifeTimeErrLow = LifeTime - (- np.mean(DriftTime) / np.log(ChargeCollectionMean - ChargeCollectionError))
    LifeTimeErrHigh = (- np.mean(DriftTime) / np.log(ChargeCollectionMean + ChargeCollectionError)) - LifeTime
    if LifeTimeErrHigh < 0: LifeTimeErrHigh = 1E10
    # print("%.2f,%.3f,%.2f,%.1f,%.3f,%.1f,%.1f" % (DriftTimeMean, ChargeCollectionMean, LifeTime, DriftTimeError, ChargeCollectionError, LifeTimeErrLow, LifeTimeErrHigh)) 

def StandardPlots(ch1, ch2):
    # ch1.Max = ch1.Max[np.where(ch1.GradTime<30)]  
    # ch2.Max = ch2.Max[np.where(ch1.GradTime<30)]  
    # ch1.TimeStamp = ch1.TimeStamp[np.where(ch1.GradTime<30)]  
    # ch1.GradTime = ch1.GradTime[np.where(ch1.GradTime<30)]  
    # ch2.GradTime = ch2.GradTime[np.where(ch1.GradTime<30)]  
    
    avg_anode = np.mean(ch1.Max)
    avg_cathode = np.mean(ch2.Max)
    # print(avg_cathode/avg_anode)


    ChargeCollection = ch1.Max / ch2.Max
    DriftTime = ch1.GradTime - ch2.GradTime
    print(" | Drift Time...", np.mean(DriftTime), np.std(DriftTime))
    print(" | Charge collection...", np.mean(ChargeCollection), np.std(ChargeCollection), np.std(ChargeCollection)/np.sqrt(len(ChargeCollection)))
    print(" | Lifetime...", -np.mean(DriftTime)/np.log(np.mean(ChargeCollection)), 'us')
    
    print(" | Drift Time...", np.median(DriftTime), np.std(DriftTime))
    print(" | Charge collection...", np.median(ChargeCollection), np.std(ChargeCollection), np.std(ChargeCollection)/np.sqrt(len(ChargeCollection)))
    print(" | Lifetime...", -np.median(DriftTime)/np.log(np.mean(ChargeCollection)), 'us')



    # print(ch1.TimeStamp[len(ch1.TimeStamp)//2], '\t', np.median(ChargeCollection), '\t', np.std(ChargeCollection)/np.sqrt(len(ChargeCollection)))
    # print(ch1.TimeStamp[len(ch1.TimeStamp)//2], '\t', np.mean(ChargeCollection), '\t', np.std(ChargeCollection)/np.sqrt(len(ChargeCollection)))



    DiffMinute = int((np.max(ch1.TimeStamp) - np.min(ch1.TimeStamp)).seconds/60.0 + 0.5)
    XTicks = int((DiffMinute/12.0 + 0.5))+1

    SavePath = '/home/fas/david_moore/aj487/purity_monitor/plots/analysis/'
    Date = datetime.datetime.now().strftime("%Y%m%d")
    Save = True 
    if Save:
        if not os.path.exists(SavePath+Date):
            os.makedirs(SavePath+Date)




    PltTime(Time=ch1.TimeStamp,Data=[ch1.Max,ch2.Max,ChargeCollection*100],Legend=['Anode','Cathode','Charge Collection [\%]'],Label='Amplitude [mV]',XTicks=XTicks,YTicks=500,YRange=[0,3000],SaveName='amp_ratio',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.GradTime,ch2.GradTime,DriftTime],Legend=['Anode','Cathode','Drift Time [$\mu$s]'],Label='Peak Time [$\mu$s]',XTicks=XTicks,YTicks=1,YRange=[0,15],SaveName='drift_time',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ch1.BaseStd,ch2.BaseStd],Legend=['Anode','Cathode'],Label='Baseline Noise [mV]',XTicks=XTicks,YTicks=2,YRange=[0,20],SaveName='baseline',Save=Save)
    # PltTime(Time=ch1.TimeStamp,Data=[ChargeCollection],Legend=[''],Label='Charge Collection',YRange=[0,2],XTicks=XTicks,YTicks=0.5,SaveName='ratio',Save=Save)
    # PltWfm(Time=ch1.Time,Data=[x+100*ii for ii,x in enumerate(ch2.Amp)],Legend=['Anode','Cathode'],XTicks=100,YTicks=200,YRange=[-50,2600],Color='k',SaveName='avg_waveform')
    # PltWfm(Time=ch2.Time,Data=[x+100*ii for ii,x in enumerate(ch2.Amp)],Legend=['Anode','Cathode'],XTicks=100,YTicks=200,YRange=[-50,2600],Color='k',SaveName='avg_waveform')
    # PltWfm(Time=ch1.Time,Data=[ch1.MeanAmp,ch2.MeanAmp],Legend=['Anode','Cathode'],XTicks=100,YTicks=50,SaveName='avg_waveform',Save=Save)
    # PltWfm(Time=ch1.Time,Data=list(ch1.AmpClean),Legend=['Anode','Cathode'],XTicks=100,YTicks=50,SaveName='avg_waveform',Color='k')
    # PltAllWfm(Time=ch2.Time,Data=list(ch1.Amp),Legend=[''],XTicks=100,YTicks=50,Save='all_waveform',Color='k')
    
    # PltWfm(Time=ch1.Time,Data=list(ch1.Amp),Legend=[''],XTicks=[50,100],YTicks=[5,10],Color='k',Save='ch1_all_waveforms')
    # PltWfm(Time=ch2.Time,Data=list(ch2.Amp),Legend=[''],XTicks=[50,100],YTicks=[5,10],Color='k',Save='ch2_all_waveforms')

    if arg.show: 
        plt.show()
    else: 
        plt.close('all')

def InitializeChannels(NumChannels, Pol=1):
    return [Wvf.WFM(ID=ii, Pol=(-1)**ii*Pol) for ii in NumChannels]

if __name__ == '__main__':
    ###### Initialize class object for each channel.
    ch1 = WFM(ID=1, Pol=1)
    ch2 = WFM(ID=2, Pol=-1)
    channels = [ch1, ch2]

    ###### Import data from HDF5 file.
    FileNames, Method = ChooseFilesToAnalyze(arg)
    for File in FileNames: 
        ImportDataFromHDF5(File, channels)

    DoAnalysis(channels)
    PrintResults(ch1, ch2)  
    StandardPlots(ch1,ch2)
    print(" | Time elapsed: ", time.process_time() , "sec")

    ###### Additional manual analysis and plotting 

import datetime
import os
import glob
import time
import sys
import WaveformClass as Wvf
import h5py
import numpy as np
import PlotFunctions as Plt
import matplotlib.pyplot as plt

def ProgressBar(it, counts):
    width = counts
    if(it==1):
        sys.stdout.write(" | [%s]" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width+1))
        sys.stdout.write("%")
    else:
        sys.stdout.write("%")
        sys.stdout.flush()
    if(it==counts):
        sys.stdout.write("\n")

def ElapsedTime():
    print(" | Time elapsed:          ", time.process_time() , "sec")

class MetaData:
    def __init__(self, InputPath):
        self.InputPath = InputPath
        # self.DateOfDataTaking = self.InputPath.split('/')[-3]
        # self.ExactTimeOfDataTaking = datetime.datetime.strptime(time.ctime(os.path.getctime(InputPath)), "%a %b %d %H:%M:%S %Y").strftime('%Y%m%d%H%M%S')
        self.DateOfProcessing = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.ProjectPath = '/project/fas/david_moore/aj487/LXeData/'
        # self.FileNames = glob.glob(self.InputPath+"*.txt")
        # self.NrOfFiles = len(self.FileNames)
        # self.OutputPath = self.CreateDirectory()
        # self.DataName = self.FileNames[0].split('/')[-2]

    def CreateDirectory(self):
        Directory = self.ProjectPath + self.DateOfDataTaking + '/'
        if not os.path.exists(Directory):
            os.makedirs(Directory)
        return Directory

class Dataset: 
    def __init__(self,  Path, ShowPlots=True, Selection='*', Pol=1, NumChannels=2):
        self.Path = Path
        self.NumChannels = NumChannels
        self.ShowPlots = ShowPlots
        self.Selection = Selection
        self.Pol = Pol 
        self.Ch = self.InitializeChannels(self.NumChannels, self.Pol)
        self.Files = glob.glob(self.Path+self.Selection)

    def RunStandardAnalysis(self): 
        # self.Files = self.ChooseFilesToAnalyze(self.Path)
        for File in self.Files: 
            self.ImportDataFromHDF5(File, self.Ch)
        self.DoAnalysis(self.Ch)
        self.ChargeCollection = self.Ch[0].Max / self.Ch[1].Max
        self.DiffMinute = int((np.max(self.Ch[0].TimeStamp) - np.min(self.Ch[0].TimeStamp)).seconds/60.0 + 0.5)
        self.XTicks = int((self.DiffMinute/12.0 + 0.5))+1
        self.Cut = np.where(self.Ch[0].BaseStd < 5)

    def InitializeChannels(self, NumChannels=2, Pol=1):
        return [Wvf.WFM(ID=ii, Pol=(-1)**ii*-1*Pol) for ii in range(1,NumChannels+1)]
    
    # def ChooseFilesToAnalyze(self, Path):
    #     return 

    def ImportDataFromHDF5(self, File, channels):
        f = h5py.File(File, 'r')  
        print(" | Filename...", File)
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

    def DoAnalysis(self, channels):
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
            ch.GetBaselineNoise(Data=ch.AmpClean)
        print(" | Time elapsed: ", time.process_time() , "sec")

    def ShowStandardPlots(self): 
        self.ShowBaselineNoise()
        self.ShowAmplitudeVsTime()
        self.ShowChargeCollection()
        self.ShowDrifttimeVsTime()

    def ShowBaselineNoise(self, Channel=-1): 
        fig = plt.figure(figsize=(8,6))
        plt.xlim(0,20)
        plt.xlabel('Baseline RMS Noise [mV]', fontsize=16)
        plt.ylabel('Counts/bin', fontsize=16)
        histmax = 0
        for ii in range(self.NumChannels):
            if ii+1 is not Channel and Channel is not -1: 
                continue
            
            h,hx,hp = plt.hist(self.Ch[ii].BaselineNoise, bins=np.arange(0.0,20.0,0.2), histtype='step', align='mid', lw=2, color=Plt.colors[ii], label=self.Ch[ii].Name)
            plt.axvline(np.median(self.Ch[ii].BaselineNoise), color=Plt.colors[ii])
            rectangle = plt.Rectangle(xy=(np.median(self.Ch[ii].BaselineNoise)-np.std(self.Ch[ii].BaselineNoise)/np.sqrt(len(self.Ch[ii].BaselineNoise)),0), 
                                    width=2*np.std(self.Ch[ii].BaselineNoise)/np.sqrt(len(self.Ch[ii].BaselineNoise)), 
                                    height=10000, 
                                    fc=Plt.colors[ii],
                                    ec=Plt.colors[ii], 
                                    alpha=0.05, 
                                    fill=True)
            plt.gca().add_patch(rectangle)
            if np.max(h) > histmax: 
                histmax = np.max(h)
        plt.ylim(0, self.RoundUpToNext(histmax, 10))
        plt.legend(loc='upper right')
        fig.tight_layout()

    def ShowAmplitudeVsTime(self, Channel=-1, YTicks=100, YMax=None): 
        if YMax is None: 
            YMax = np.max([np.max(self.Ch[ii].Max) for ii in range(self.NumChannels)])
#         YMax = self.RoundUpToNext(YMax, 100)
#         print('Max: ', YMax)
#         YTicks = self.RoundDownToNext(YMax/5, 10)
        print('Ticks: ', YTicks)
        Plt.PltTime(Time=self.Ch[0].TimeStamp[self.Cut],
                    Data=[self.Ch[0].Max[self.Cut], self.Ch[1].Max[self.Cut], self.ChargeCollection[self.Cut]*100],
                    Legend=['Anode','Cathode','Charge Collection [\%]'],
                    Label='Amplitude [mV]',
                    XTicks=self.XTicks,
                    YTicks=YTicks,
                    YRange=[0,YMax],
                    SaveName='amp_ratio',
                    Title='',
                    Save=False)

    def ShowDrifttimeVsTime(self, Channel=-1): 
        YMax = np.max([np.max(self.Ch[ii].GradTime) for ii in range(self.NumChannels)])
        YMax = self.RoundUpToNext(YMax, 10)
        YTicks = self.RoundDownToNext(YMax/5, 1)
        Plt.PltTime(Time=self.Ch[0].TimeStamp,
                    Data=[self.Ch[0].GradTime, self.Ch[1].GradTime, self.Ch[0].GradTime - self.Ch[1].GradTime],
                    Legend=['Anode','Cathode','Drift Time'],
                    Label='Time since Trigger [$\mu$s]',
                    XTicks=self.XTicks,
                    YTicks=YTicks,
                    YRange=[0,YMax],
                    SaveName='drift_time',
                    Save=False)

    def ShowChargeCollection(self, Channel=-1): 
        Mean = np.mean(self.ChargeCollection)
        Median = np.median(self.ChargeCollection)
        Err = np.std(self.ChargeCollection)/np.sqrt(len(self.ChargeCollection))       

        XMin = self.RoundDownToNext(np.min(self.ChargeCollection), 1)
        XMax = self.RoundUpToNext(np.max(self.ChargeCollection), 1)
        fig = plt.figure(figsize=(8,6))
        
        h,hx,hp = plt.hist(self.ChargeCollection, bins=np.arange(XMin,XMax,0.01), histtype='step', align='mid', lw=2)

        plt.xlim(XMin, self.RoundUpToNext(np.max(self.ChargeCollection), 1))
        plt.ylim(0, self.RoundUpToNext(np.max(h), 10))

        plt.xlabel('Charge Collection', fontsize=16)
        plt.ylabel('Counts/bin', fontsize=16)
        plt.axvline(Mean, color='red', )
        plt.title('Charge Collection = %.4f ± %.4f' % (Mean, Err))
        rectangle = plt.Rectangle((Mean-Err,0), 2*Err, 10000, fc='red',ec="red", alpha=0.2, fill=True,)
        plt.gca().add_patch(rectangle)
        fig.tight_layout()

    def RoundUpToNext(self, Num, Ceil): 
        return int(np.ceil(Num / float(Ceil))) * float(Ceil)

    def RoundDownToNext(self, Num, Floor): 
        return int(np.floor(Num / float(Floor))) * float(Floor)

import datetime
import glob
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt

import PlotFunctions as Plt
import Waveform as Wvf

class Dataset: 
    def __init__(self,  Path, ShowPlots=True, Selection='*', Pol=1, NumChannels=2):
        self.Path = Path
        self.NumChannels = NumChannels
        self.ShowPlots = ShowPlots
        self.Selection = Selection
        self.Pol = Pol 
        self.Ch = self.InitializeChannels(self.NumChannels, self.Pol)
        self.Files = glob.glob(self.Path+self.Selection)

    def RunStandardAnalysis(self, NoiseDataset=None): 
        # self.Files = self.ChooseFilesToAnalyze(self.Path)
        for File in self.Files: 
            self.ImportDataFromHDF5(File, self.Ch)
        self.DoAnalysis(self.Ch, NoiseDataset=NoiseDataset)
        self.ChargeCollection = self.Ch[0].Max / self.Ch[1].Max
        self.DiffMinute = int((np.max(self.Ch[0].TimeStamp) - np.min(self.Ch[0].TimeStamp)).seconds/60.0 + 0.5)
        self.XTicks = int((self.DiffMinute/12.0 + 0.5))+1
        self.NoiseCut = 1000
        self.Cut = np.where(self.Ch[0].BaseStd < self.NoiseCut)
        self.InverseCut = np.where(self.Ch[0].BaseStd > self.NoiseCut)


    def InitializeChannels(self, NumChannels=2, Pol=1):
        return [Wvf.Waveform(ID=ii, Pol=(-1)**ii*-1*Pol) for ii in range(1,NumChannels+1)]
    
    # def ChooseFilesToAnalyze(self, Path):
    #     return 

    def ImportDataFromHDF5(self, File, channels):
        f = h5py.File(File, 'r')  
        # print(" | Filename...", File)
        Keys = list(f.keys())
        for ch in channels:
            ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
            ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
            Group = f.get(ch.ChName)
            GroupKeys = Group.keys()
            ch.Files.append(len(GroupKeys))
            # print(" | Number of files in ch%d...\t %d/%d" % (ch.ID, ch.Files[-1], np.sum(ch.Files)))
            for key in GroupKeys:
                ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale * ch.Pol)
                ch.TimeStamp.append(datetime.datetime.strptime((f.attrs['Date']+Group.get(key).attrs["TimeStamp"]), '%Y%m%d%H%M%S'))
        f.close()
            

    def DoAnalysis(self, channels, NoiseDataset=None):
    ###### Basic analysis: baseline subtraction, waveform averaging, obtaining fourier spectra, frequency bandpass filter and finding extrema.
        Print = False 
        for ii, ch in enumerate(channels):
            # print(" | Processing data in channel %d..." % (ch.ID))
            ch.GetSampling()
            ch.Amp = [x for _, x in sorted(zip(ch.TimeStamp, ch.Amp))]
            ch.Amp = np.array(ch.Amp)
            ch.TimeStamp = np.array(sorted(ch.TimeStamp))

            # ch.TimeStamp = np.array(ch.TimeStamp)
            ch.Amp = ch.SubtractBaseline(Data=ch.Amp, state=Print)
            ch.Amp = ch.RemoveNoise(Data=ch.Amp, HighPass=80000, state=Print)
            if NoiseDataset is not None:
                for jj,amp in enumerate(ch.Amp):
                    ch.Amp[jj] =  ch.Amp[jj]-np.mean(NoiseDataset.Ch[ii].Amp,axis=0)

            # ch.RunFit(Data=ch.Amp)
            ch.GetAllMaxima(Data=ch.Amp, state=Print)
            # ch.FindMaxGradient(Data=ch.Amp ,state=Print)
            ch.GetDriftTime(Data=ch.Amp, Threshold=0.1)
            ch.GetIntegral(Data=ch.Amp, state=Print)
            # ch.GetBaselineNoise(Data=ch.Amp)

        # print(" | Time elapsed: ", time.process_time() , "sec")

    def ShowStandardPlots(self): 
        self.ShowBaselineNoise()
        self.ShowAmplitudeVsTime()
        self.ShowChargeCollection()
        self.ShowDrifttimeVsTime()

    def ShowBaselineNoise(self, Channel=-1, BinMax=20): 
        fig = plt.figure(figsize=(8,6))
        plt.xlim(0,BinMax)
        plt.xlabel('Baseline RMS Noise [mV]')
        plt.ylabel('Counts/bin')
        histmax = 0
        for ii in range(self.NumChannels):
            if ii+1 is not Channel and Channel != -1: 
                continue
            
            h,hx,hp = plt.hist(self.Ch[ii].BaseStd[self.Cut], bins=np.arange(0.0,BinMax,0.2), histtype='step', align='mid', lw=2, color=Plt.colors[ii], label=self.Ch[ii].Name)
            plt.axvline(np.median(self.Ch[ii].BaseStd[self.Cut]), color=Plt.colors[ii])
            rectangle = plt.Rectangle(xy=(np.median(self.Ch[ii].BaseStd[self.Cut])-np.std(self.Ch[ii].BaseStd[self.Cut])/np.sqrt(len(self.Ch[ii].BaseStd[self.Cut])),0), 
                                    width=2*np.std(self.Ch[ii].BaseStd[self.Cut])/np.sqrt(len(self.Ch[ii].BaseStd[self.Cut])), 
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
        Plt.PltChargeVsTime(Time=self.Ch[0].TimeStamp[self.Cut],
                    Data=[self.Ch[0].Max[self.Cut], self.Ch[1].Max[self.Cut]],
                    CC=self.ChargeCollection[self.Cut]*100,
                    Legend=['Anode','Cathode','Charge Collection'],
                    Label='Amplitude [mV]',
                    XTicks=self.XTicks,
                    YTicks=YTicks,
                    YRange=[0,YMax],
                    Title='')

    def ShowIntegralVsTime(self, Channel=-1, YTicks=100, YMax=None): 
        
        self.Ch[0].Integral /= np.max(self.Ch[0].Integral)
        self.Ch[1].Integral /= np.max(self.Ch[1].Integral)
        if YMax is None: 
            YMax = np.max([np.max(self.Ch[ii].Integral) for ii in range(self.NumChannels)])
        Plt.PltChargeVsTime(Time=self.Ch[0].TimeStamp[self.Cut],
                    Data=[self.Ch[0].Integral[self.Cut], self.Ch[1].Integral[self.Cut]],
                    CC=self.Ch[0].Integral[self.Cut]/self.Ch[1].Integral[self.Cut]*100,
                    Legend=['Anode','Cathode','Charge Collection'],
                    Label='Integral [a.u.]',
                    XTicks=self.XTicks,
                    YTicks=YTicks,
                    YRange=[0,YMax],
                    Title='')

    def ShowDrifttimeVsTime(self, YMax=0, Channel=-1): 
        # YMax = np.max([np.max(self.Ch[ii].DriftTime) for ii in range(self.NumChannels)])
        YMax = self.RoundUpToNext(YMax, 10)
        YTicks = self.RoundDownToNext(YMax/5, 1)
        self.DriftTime = self.Ch[0].DriftTime - self.Ch[1].DriftTime
        DriftCut = np.where(self.DriftTime > 0)
        Plt.PltTime(Time=self.Ch[0].TimeStamp[DriftCut],
                    Data=[self.Ch[0].DriftTime[DriftCut], self.Ch[1].DriftTime[DriftCut], self.Ch[0].DriftTime[DriftCut] - self.Ch[1].DriftTime[DriftCut]],
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

        plt.xlabel('Charge Collection')
        plt.ylabel('Counts/bin')
        plt.axvline(Mean, color='red', )
        plt.title('Charge Collection = %.4f ± %.4f' % (Mean, Err))
        rectangle = plt.Rectangle((Mean-Err,0), 2*Err, 10000, fc='red',ec="red", alpha=0.2, fill=True,)
        plt.gca().add_patch(rectangle)
        fig.tight_layout()

    def RoundUpToNext(self, Num, Ceil): 
        return int(np.ceil(Num / float(Ceil))) * float(Ceil)

    def RoundDownToNext(self, Num, Floor): 
        return int(np.floor(Num / float(Floor))) * float(Floor)

import numpy as np
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from scipy.optimize import curve_fit

from PlotFunctions import *
from FitFunctions import * 
from MatchedFilter import MatchedFilter 

class WFM:
    def __init__(self, ID, Pol, VScale = "m", TScale = "u"):
        self.NameDict = {1:'Anode', -1:'Cathode'}
        self.TimeStamp = []
        self.ID = ID
        self.Pol = Pol
        self.ChName = 'ch%d' % self.ID
        self.Name = self.NameDict[self.Pol]
        self.Files = []
        self.Samples = 0
        self.Sampling = 0
        self.BaseCounts = 2500
        self.Amp = []
        self.AmpClean = []
        self.AmpSubtract = []
        self.Baseline = []
        self.BaseStd = []
        self.MeanAmp = []
        self.Time = []
        self.Trigger = []
        self.Max = []
        self.MaxT = []
        self.AmpFFT = []
        self.TimeFFT = []
        self.VScale = self.Scale(units = VScale)
        self.TScale = self.Scale(units = TScale)
        self.Plot = False
        self.Integral = []
        self.Voltages = []

    def Scale(self, units):
        if(units=="n"): return 1000000000.0
        if(units=="u"): return 1000000.0
        if(units=="m"): return 1000.0
        if(units=="1"): return 1.0

    def FindTimeBin(self, bin):
        return np.abs((np.asarray(self.Time)-bin)).argmin()

    def SetPolarity(self):
        if self.Label == 'Cathode':
            self.Pol = -1
            print(" | This is a cathode signal...")
        if self.Label == 'Cathode Grid':
            self.Pol = 1
            print(" | This is not a cathode signal...")
        if self.Label == 'Anode':
            self.Pol = +1
            print(" | This is an anode signal...")

    def GetSampling(self, state=False):
        if(state): print(" | Get sampling...")
        self.Sampling = self.TScale/abs(self.Time[0]-self.Time[1])
        self.Samples = len(self.Time)

    def SubtractBaseline(self, state=False):
        if(state): print(" | Subtracting baseline...")
        self.BaseCounts = self.FindTimeBin(-50)
        self.BaseStd = []
        self.Baseline = []
        for i in range(np.sum(self.Files)):
            self.BaseStd.append(np.std(self.AmpClean[i][self.FindTimeBin(-900):self.FindTimeBin(-10)]))
            self.Baseline.append(np.average(self.AmpClean[i][self.FindTimeBin(-50):self.FindTimeBin(0)]))
            # self.Baseline.append(np.average(self.Amp[i][:self.BaseCounts]))
            self.AmpClean[i] = self.AmpClean[i] - self.Baseline[i]
        self.BaseStd = np.array(self.BaseStd)
        self.Baseline = np.array(self.Baseline)

    def ApplyCut(self, Cut, state=False):
        Cut = Cut[0]
        if(state): print(" | Applying cuts... %d out %d left" % (len(Cut),len(self.Amp)))
        self.Amp = self.Amp[Cut]
        self.AmpClean = self.AmpClean[Cut]
        self.Max = self.Max[Cut]
        self.MaxT = self.MaxT[Cut]
        self.GradTime = self.GradTime[Cut]
        self.BaseStd = self.BaseStd[Cut]
        self.TimeStamp = self.TimeStamp[Cut]

    def LinearFunction(self,x,slope,yoffset):
        return x*slope + yoffset

    def SubtractLinearBackground(self, Time, Data, state=False): 
        if(state): print(" | Subtracting linear background...")
        for data in Data:
            UpLimit = -100
            init_vals = [0,0]
            best_vals, covar = curve_fit(self.LinearFunction, Time[np.where(Time<UpLimit)], data[np.where(Time<UpLimit)], p0=init_vals)
            self.AmpSubtract.append(np.array(data - self.LinearFunction(Time, *best_vals)))
        self.AmpSubtract = np.array(self.AmpSubtract)

    def SubtractFunction(self, data, fit, start, end, state=False):
        if(state): print(" | Subtracting baseline...")
        new = np.subtract(data,fit, where=((self.Time > start) & (self.Time < end)))
        return new

    def FindMaxGradient(self,Data,state=False):
        if(state): print(" | Calculating gradient along waveform...")
        Spacing = 50
        Bins = np.linspace(0,len(self.Time[::Spacing])*Spacing,len(self.Time[::Spacing]))
        self.GradTime = []
        for data in Data: 
            Gradient = np.gradient(data[::Spacing])
            MaxGradient = np.where(Gradient==np.max(Gradient[len(Gradient)//2:len(Gradient)*3//5]))
            MaxGradientBin = int(Bins[MaxGradient])
            GradTime = self.Time[MaxGradientBin]
            if GradTime < 0:
                print('STOPPP')
                PltWfm(Time=self.Time,Data=[data],Legend=[''],XTicks=100,YTicks=50,SaveName='test',Save=False)
                plt.show()
            self.GradTime.append(self.Time[MaxGradientBin])
        self.GradTime = np.array(self.GradTime)

    def GetAverageWfm(self, Data, state=False):
        if(state): print(" | Getting average waveform...")
        self.MeanAmp = np.mean(Data, axis=0)

    def GetIntegral(self, Data, state=False):
        if(state): print(" | Getting average waveform...")
        for i in range(np.sum(self.Files)):
            self.Integral.append(np.sum(Data[i][self.FindTimeBin(0):self.FindTimeBin(200)])/100000.0)
        self.Integral = np.array(self.Integral)
        
    def GetAllMaxima(self, data, state=False):
        self.Max = []
        self.MaxT = []
        if(state): print(" | Getting extrema of individual files...")
        for i in range(np.sum(self.Files)):
            self.Max.append(np.max(data[i,self.FindTimeBin(0):self.FindTimeBin(150)]))
            self.MaxT.append(self.Time[np.where(data[i]==self.Max[i])[0][0]])
        self.Max = np.array(self.Max)
        self.MaxT = np.array(self.MaxT)

    def GetAllFFT(self, state=False):
        if(state): print(" | Getting Fourier spectra...")
        T = 1.0 / self.Sampling
        self.TimeFFT = np.linspace(0.0, 1.0/(2.0*T), self.Samples/2)[1:self.Samples//2]
        for i in range(np.sum(self.Files)):
            self.AmpFFT.append(np.abs(fft(self.Amp[i])[1:self.Samples//2]).tolist())

    def RemoveNoise(self,LowCut, HighCut, Order, state=False):
        self.AmpClean = [] 
        if(state): print(" | Removing noise...")
        for i in range(np.sum(self.Files)):
            self.AmpClean.append(self.butter_bandpass_filter(self.Amp[i], LowCut, HighCut, self.Sampling, Order))
        self.AmpClean = np.array(self.AmpClean)

    def butter_bandpass(self, lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], btype='bandpass', output='sos')
        return sos

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
        sos = self.butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y

    def RemoveNoiseSingle(self, data, LowCut, HighCut, Order):
        return self.butter_bandpass_filter(data, LowCut, HighCut, self.Sampling, Order).tolist()

    def ShapeGaussian(self, data, sigma):
        gaussian = []
        for x in self.Time:
            gaussian.append(np.exp(-((x/sigma)**2)/2)*(1/np.sqrt(2*np.pi*sigma*sigma)))
        result = np.convolve(data, gaussian, mode="same")
        return result.tolist()

    def GetBaselineNoise(self, Data):
        self.BaselineNoise = [np.sqrt(np.mean(data[self.FindTimeBin(-900):self.FindTimeBin(-100)]**2)) for data in Data] 

    def GetDriftTime(self, Data, Threshold=0.1): 
        self.DriftTime = []
        for ii,data in enumerate(Data): 
            if self.Pol == 1:
                ThresholdVal = (1-Threshold+0.05) * self.Max[ii]
            elif self.Pol == -1: 
                ThresholdVal = Threshold * self.Max[ii]
            try:
                BinsOverThreshold = np.where(data[self.FindTimeBin(0):self.FindTimeBin(self.MaxT[ii])] > ThresholdVal )[0][0]
            except IndexError as error:
                print(data[self.FindTimeBin(7):self.FindTimeBin(self.MaxT[ii])], ThresholdVal)
            self.DriftTime.append(self.Time[self.FindTimeBin(0)+BinsOverThreshold])
        self.DriftTime = np.array(self.DriftTime)
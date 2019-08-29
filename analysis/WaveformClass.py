import numpy as np
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from scipy.optimize import curve_fit

from PlotFunctions import *
from FitFunctions import * 
from MatchedFilter import MatchedFilter 

class WFM:
    def __init__(self, ID, VScale = "m", TScale = "u"):
        self.TimeStamp = []
        self.ID = ID
        self.ChName = 'ch%d' % self.ID
        self.Files = 0
        self.counter = 0
        self.Pol = 1
        self.Samples = 0
        self.Sampling = 0
        self.BaseCounts = 2500
        self.Amp = [[] for x in range(self.Files)]
        self.AmpClean = [[] for x in range(self.Files)]
        self.AmpClean = [[] for x in range(self.Files)]
        self.AmpSubtract = [[] for x in range(self.Files)]
        self.Baseline = []
        self.BaseStd = []
        self.MeanAmp = []
        self.Peak = []
        self.PeakTime = []
        self.Time = []
        self.Trigger = []
        self.Max = []
        self.MaxT = []
        self.AmpFFT = []
        self.TimeFFT = []
        self.VScale = self.Scale(units = VScale)
        self.TScale = self.Scale(units = TScale)
        self.Plot = False
        self.CreateTime = []
        self.Integral = []

    def count(self):
        self.counter = self.counter + 1

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
        for i in range(self.Files):
            self.BaseStd.append(np.std(self.Amp[i][:self.BaseCounts]))
            self.Baseline.append(np.average(self.Amp[i][self.FindTimeBin(-50):self.FindTimeBin(0)]))
            # self.Baseline.append(np.average(self.Amp[i][:self.BaseCounts]))
            self.Amp[i] = self.Amp[i] - self.Baseline[i]
        self.BaseStd = np.array(self.BaseStd)
        self.Baseline = np.array(self.Baseline)

    def LinearFunction(self,x,slope,yoffset):
        return x*slope + yoffset

    def SubtractLinearBackground(self, Time, Data, state=False): 
        if(state): print(" | Subtracting linear background...")
        for data in Data:
            UpLimit = -100
            init_vals = [0,0]
            best_vals, covar = curve_fit(self.LinearFunction, Time[np.where(Time<UpLimit)], data[np.where(Time<UpLimit)], p0=init_vals)
            self.AmpSubtract.append(np.array(data - self.LinearFunction(Time, best_vals[0], best_vals[1])))
            # plt.xlim(-1000,1000)
            # plt.scatter(Time,Data,s=0.1, c='black')
            # plt.plot(Time,Data-self.LinearFunction(Time, best_vals[0], best_vals[1]),linewidth=0.02, c='yellow')
            # plt.plot(Time, self.LinearFunction(Time, best_vals[0], best_vals[1]), label=r'f(x)=$%f \cdot x + %f$' % (best_vals[0], best_vals[1]), linewidth=1.5, c='red')
            # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
            # plt.show()
        self.AmpSubtract = np.array(self.AmpSubtract)

    def SubtractFunction(self, data, fit, start, end, state=False):
        if(state): print(" | Subtracting baseline...")
        new = np.subtract(data,fit, where=((self.Time > start) & (self.Time < end)))
        return new

    def GetAverageWfm(self, Data, state=False):
        if(state): print(" | Getting average waveform...")
        self.MeanAmp = np.mean(Data, axis=0)

    def GetIntegral(self, Data, state=False):
        if(state): print(" | Getting average waveform...")
        for i in range(self.Files):
            self.Integral.append(np.sum(Data[i][self.FindTimeBin(0):]))
        self.Integral = np.array(self.Integral)
        
    def GetAllMaxima(self, data, state=False):
        self.Max = []
        self.MaxT = []
        if(state): print(" | Getting extrema of individual files...")
        if(self.Pol==1):
            for i in range(self.Files):
                self.Max.append(np.max(data[i,self.FindTimeBin(20):]))
                self.MaxT.append(self.Time[np.where(data[i]==self.Max[i])[0][0]])
        else:
            for i in range(self.Files):
                self.Max.append(np.min(data[i]))
                self.MaxT.append(self.Time[np.where(data[i]==self.Max[i])[0][0]])
        self.Max = np.array(self.Max)
        self.MaxT = np.array(self.MaxT)

    def GetAllFFT(self, state=False):
        if(state): print(" | Getting Fourier spectra...")
        T = 1.0 / self.Sampling
        self.TimeFFT = np.linspace(0.0, 1.0/(2.0*T), self.Samples/2)[1:self.Samples//2]
        for i in range(self.Files):
            self.AmpFFT.append(np.abs(fft(self.Amp[i])[1:self.Samples//2]).tolist())

    def RemoveNoise(self,LowCut, HighCut, Order, state=False):
        self.AmpClean = [] 
        if(state): print(" | Removing noise...")
        for i in range(self.Files):
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

    def GetPeak(self, data):
        if(self.Pol==1):
            self.Peak = max(data)
        else:
            self.Peak = min(data)
        self.PeakTime = self.Time[data.index(self.Peak)]
        return self.Peak, self.PeakTime







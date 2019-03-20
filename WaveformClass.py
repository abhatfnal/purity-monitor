import sys
import numpy as np
import scipy.special as spl
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from scipy.ndimage.filters import gaussian_filter
from plots import PPltWfm, PltWfm, plot_single_fft

import ROOT
# ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptTitle(111);
fac = ROOT.TMath.Factorial

class WFM:
    def __init__(self, Directory, ID, VScale = "m", TScale = "u"):
        self.TimeStamp = []
        self.Directory = Directory
        self.ID = ID
        self.ChName = 'ch%d' % self.ID
        self.Files = 0
        self.counter = 0
        self.Pol = 1
        self.Samples = 0
        self.Sampling = 0
        self.BaseCounts = 2500
        self.Amp = [[] for x in range(self.Files)]
        self.Baseline = []
        self.BaseStd = []
        self.MeanAmp = []
        self.Peak = []
        self.PeakTime = []
        self.Time = []
        self.Trigger = []
        self.Max = []
        self.Int = []
        self.MaxT = []
        self.AmpFFT = []
        self.AmpClean = []
        self.TimeFFT = []
        self.VScale = self.Scale(units = VScale)
        self.TScale = self.Scale(units = TScale)
        self.Label = self.GetChannelName()
        self.Plot = False
        self.CreateTime = []

    def count(self):
        self.counter = self.counter + 1

    def Scale(self, units):
        if(units=="n"): return 1000000000.0
        if(units=="u"): return 1000000.0
        if(units=="m"): return 1000.0
        if(units=="1"): return 1.0

    def GetChannelName(self):
        label = ''
        if(self.ID == 2):
            if 'cathode' in self.Directory:
                label = 'Cathode'
            if 'anode' in self.Directory:
                label = 'Anode'
            if 'anodegrid' in self.Directory:
                label = 'Anode Grid'
            if 'cathodegrid' in self.Directory:
                label = 'Cathode Grid'
            if 'noise_' in self.Directory:
                label = 'Noise'
        else:
            label = 'Trigger'
        return label

    def FindTimeBin(self, bin):
        return np.abs((np.asarray(self.Time)-bin)).argmin()

    def SetPolarity(self):
        if self.Label == 'Cathode':
            self.Pol = -1
            print " | This is a cathode signal..."
        if self.Label == 'Cathode Grid':
            self.Pol = 1
            print " | This is not a cathode signal..."
        if self.Label == 'Anode':
            self.Pol = +1
            print " | This is an anode signal..."

    def GetSampling(self, state=False):
        if(state): print " | Get sampling..."
        self.Sampling = self.TScale/abs(self.Time[0]-self.Time[1])
        self.Samples = len(self.Time)

    def SubtractBaseline(self, state=False):
        if(state): print " | Subtracting baseline..."
        self.BaseCounts = self.FindTimeBin(-50)
        for i in range(self.Files):
            self.BaseStd.append(np.std(self.Amp[i][:self.BaseCounts]))
            self.Baseline.append(np.average(self.Amp[i][:self.BaseCounts]))
            self.Amp[i] = [(x - self.Baseline[i]) for x in self.Amp[i]]

    def SubtractFunction(self, data, fit, start, end, state=False):
        if(state): print " | Subtracting baseline..."
        DataMinusFit = [(x-y) for x,y in zip(data[start:end],fit[start:end])]
        new = data[:start]
        new.extend(DataMinusFit)
        new.extend(data[end:])
        return new

    def GetAverageWfm(self, state=False):
        if(state): print " | Getting average waveform..."
        self.MeanAmp = np.mean(self.Amp, axis=0)

    def GetAllMaxima(self, data, state=False):
        if(state): print " | Getting extrema of individual files..."
        for i in range(self.Files):
            if(self.Pol==1):
                self.Max.append(max(data[i]))
            else:
                self.Max.append(min(data[i]))
            self.MaxT.append(self.Time[data[i].index(self.Max[i])])

    def GetAllFFT(self, state=False):
        if(state): print " | Getting Fourier spectra..."
        T = 1.0 / self.Sampling
        self.TimeFFT = np.linspace(0.0, 1.0/(2.0*T), self.Samples/2)[1:self.Samples//2]
        for i in range(self.Files):
            self.AmpFFT.append(np.abs(fft(self.Amp[i])[1:self.Samples//2]).tolist())

    def RemoveNoise(self,LowCut, HighCut, Order, state=False):
        if(state): print " | Removing noise..."
        for i in range(self.Files):
            self.AmpClean.append(self.butter_bandpass_filter(self.Amp[i], LowCut, HighCut, self.Sampling, Order).tolist())

    def butter_bandpass(self, lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

    def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
        sos = self.butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y

    def RemoveNoiseSingle(self, data, LowCut, HighCut, Order):
        return self.butter_bandpass_filter(data, LowCut, HighCut, self.Sampling, Order).tolist()

    def FitExponential(self, data, start, end, repeats, state=False):
        if(state): print " | Fitting exponential to decaying edge..."
        hist = ROOT.TH1F("waveform", "waveform", self.Samples, self.Time[0],self.Time[-1])
        for i in range(len(data)):
            hist.AddBinContent(i, data[i])
        f = ROOT.TF1("f","[0]*exp(-[1]*(x-[2]))+[3]", start, end)
        f.SetLineColor(ROOT.kRed);
        f.SetNpx(100000);

        ParNames = ['Amplitude', 'Fall Time', 'X Offset', 'Y Offset']
        ParValues = [1.0, 0.001, 0.0, 0.0]
        ParLow = [0, 0, -10, -10]
        ParHigh = [10000, 10000, 10, 10]

        for i in range(len(ParNames)):
            f.SetParameter(i, ParValues[i])
            f.SetParLimits(i, ParLow[i], ParHigh[i])
            f.SetParName(i, ParNames[i])

        for i in range(repeats):
            hist.Fit("f","RME","");
        fit = self.exponential_func(np.asarray(self.Time), f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3))
        return fit

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

    def GetFitAmp(self, data, num):
        ScaleUp = 1.5
        ScaleLow = 0.5
        if(self.Pol == -1):
            amp = np.min(data)
            low = np.min(data)*ScaleUp
            high = np.min(data)*ScaleLow
        else:
            amp = np.max(data)
            low = np.max(data)*ScaleLow
            high = np.max(data)*ScaleUp
        if(num==0):
            return amp
        if(num==1):
            return low
        if(num==2):
            return high

    def FitFullCurve(self, data, start, end, repeats, state=False):
        if(state): print " | Fitting function to total curve..."
        hist = ROOT.TH1F("waveform", "waveform", self.Samples, self.Time[0],self.Time[-1])
        for i in range(len(data)):
            hist.AddBinContent(i, data[i])
        f = ROOT.TF1("f", "([3]/2)*exp(0.5*([1]/[2])^2 -(x-[4])/[2])*erfc(([1]/[2]-(x-[4])/[1])/sqrt(2))+[0]", start, end)
        f.SetLineColor(ROOT.kRed)
        f.SetNpx(10000)
    	f.SetNumberFitPoints(10000);

        ParNames = ['Baseline', 'Rise Time', 'Fall Time', 'Amplitude', 'Peak Time']
        ParValues = [0.0, 1.0, 100.0, self.GetFitAmp(data,0), 10.0]
        ParLow = [-1, 0, 100.0, self.GetFitAmp(data,1), 0]
        ParHigh = [1, 15, 100000.0, self.GetFitAmp(data,2), 1000]

        for i in range(len(ParNames)):
            f.SetParameter(i, ParValues[i])
            f.SetParLimits(i, ParLow[i], ParHigh[i])
            f.SetParName(i, ParNames[i])

        for i in range(repeats-1):
            print (" | Fit repition... %d" % i)
            hist.Fit("f", "REQM", "")
        hist.Fit("f", "REM", "")

        FitParameters = []
        print ' | Extremum:', f.GetMaximum(), f.GetMinimum()
        print ' | Position:', f.GetMaximumX(), f.GetMinimumX()
        print " | Reduced chi square...", f.GetChisquare()/f.GetNDF()
        for i in range(len(ParNames)):
            print " | Fit parameters...", ParNames[i], f.GetParameter(i)
            FitParameters.append(f.GetParameter(i))

        print ' | Graph printed...Press any key to continue...'
        raw_input()
        return FitParameters

    def FitFullCurveDouble(self, data, start, end, repeats, state=False):
        if(state): print " | Fitting function to total curve..."
        hist = ROOT.TH1F("waveform2", "waveform2", self.Samples, self.Time[0],self.Time[-1])
        for i in range(len(data)):
            hist.AddBinContent(i, data[i])

        f = ROOT.TF1("f2", "[3]/2*exp(0.5*([1]/[2])^2 - (x-[4])/[2]) * erfc(([1]/[2]-(x-[4])/[1])/sqrt(2))+[0] + [8]/2*exp(0.5*([6]/[7])^2 - (x-[9])/[7]) * erfc(([6]/[7]-(x-[9])/[6])/sqrt(2))+[5]", start, end)
        f.SetLineColor(ROOT.kRed)
        f.SetNpx(100000)

        ParNames = ['Baseline', 'Rise Time', 'Fall Time', 'Amplitude', 'Peak Time', 'Baseline 2', 'Rise Time 2', 'Fall Time 2', 'Amplitude 2', 'Peak Time 2']
        ParValues = [0.0, 1.0, 100.0, self.GetFitAmp(data,0), 10.0, 0.0, 1.0, 100.0, self.GetFitAmp(data,0), 10.0]
        ParLow = [-1, 0, 100.0, self.GetFitAmp(data,1), 0, -1, 0, 100.0, self.GetFitAmp(data,1), 0]
        ParHigh = [1, 15, 100000.0, self.GetFitAmp(data,2), 1000, 1, 15, 100000.0, self.GetFitAmp(data,2), 1000]

        ParValues[:5] = self.FitFullCurve(data, -100, 70, 5)
        ParValues[5:10] = self.FitFullCurve(data, 90, 300, 5)
        ParLow = 0.5*np.asarray(ParValues)
        ParHigh = 1.5*np.asarray(ParValues)
        print len(ParLow)
        print len(ParHigh)
        print len(ParValues)
        print ParValues
        print 'seperate fits done'

        for i in range(len(ParNames)):
            print i
            f.SetParameter(i, ParValues[i])
            f.SetParLimits(i, ParLow[i], ParHigh[i])
            f.SetParName(i, ParNames[i])

        for i in range(repeats-1):
            print (" | Fit repition... %d" % i)
            hist.Fit("f2", "REMQ", "")
        hist.Fit("f2", "REM", "")

        FitParameters = []
        print ' | Extremum:', f.GetMaximum(), f.GetMinimum()
        print ' | Position:', f.GetMaximumX(), f.GetMinimumX()
        print " | Reduced chi square...", f.GetChisquare()/f.GetNDF()
        for i in range(len(ParNames)):
            print " | Fit parameters...", ParNames[i], f.GetParameter(i)
            FitParameters.append(f.GetParameter(i))

        print ' | Graph printed...Press any key to continue...'
        raw_input()
        return FitParameters

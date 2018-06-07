import numpy as np
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, sosfilt
from scipy.ndimage.filters import gaussian_filter
from plots import PPltWfm, plot_single_waveform, plot_single_fft

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptTitle(111);
fac = ROOT.TMath.Factorial

class WFM:
    def __init__(self, Files, Pol, Label, VScale = "m", TScale = "u"):
        self.counter = 0
        self.Files = Files
        self.Pol = Pol
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
        self.Max = []
        self.Int = []
        self.MaxT = []
        self.AmpFFT = []
        self.AmpClean = []
        self.TimeFFT = []
        self.VScale = self.Scale(units = VScale)
        self.TScale = self.Scale(units = TScale)
        self.label = Label

    def count(self):
        self.counter = self.counter + 1

    def Scale(self, units):
        if(units=="n"): return 1000000000.0
        if(units=="u"): return 1000000.0
        if(units=="m"): return 1000.0
        if(units=="1"): return 1.0

    def GetSampling(self):
        self.Sampling = self.TScale/abs(self.Time[0]-self.Time[1])
        self.Samples = len(self.Time)

    def SubtractBaseline(self, state=False):
        if(state): print " | Subtracting baseline..."
        for i in range(self.Files):
            self.BaseStd.append(np.std(self.Amp[i][:self.BaseCounts]))
            self.Baseline.append(np.average(self.Amp[i][:self.BaseCounts]))
            self.Amp[i] = [(x - self.Baseline[i]) for x in self.Amp[i]]

    def GetAverageWfm(self, state=False):
        if(state): print " | Getting average waveform..."
        for i in range(self.Samples):
            self.MeanAmp.append(np.average([self.Amp[j][i] for j in range(len(self.Amp))]))
        _,_ = self.GetPeak(self.MeanAmp, self.Pol)

    def GetAllMaxima(self, data, state=False):
        if(state): print " | Getting extrema of individual files..."
        for i in range(self.Files):
            if(self.Pol==1):
                self.Max.append(max(data[i]))
            else:
                self.Max.append(min(data[i]))
            self.MaxT.append(self.Time[data[i].index(self.Max[i])])
            # print " | ", i, self.MaxT[i], self.Max[i]

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

        f.SetParameter(0,1.0);
        f.SetParLimits(0,-1000,1000);
        f.SetParName(0, "Amp");
        f.SetParameter(1,.001);
        f.SetParLimits(1,-10,10);
        f.SetParName(1, "Tau");
        f.SetParameter(2,0);
        f.SetParLimits(2,-100,100);
        f.SetParName(2, "XOff");
        f.SetParameter(3,0);
        f.SetParLimits(3,-100,100);
        f.SetParName(3, "YOff");

        for i in range(repeats):
            hist.Fit("f","RMEQ","");
        fit = self.exponential_func(np.asarray(self.Time), f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3))
        # PPltWfm(self.Time, data, fit, 'Cathode', 'Fit','Time [$\mu$s]', 'Amplitude [mV]',scale=1.2,xlim=self.Time[0],xlim2=self.Time[-1],ylim=min(data)*1.2,ylim2=max(data)*1.2)
        return start, fit[self.Time.index(start)]

    def exponential_func(self, x, a, b, c, d):
        return a*np.exp(-b*(x-c))+d

    def double_exponential_func(self, x, a, b,  d, e, g):
        return a*np.exp(-b*(x))+d*np.exp(-e*(x))+g

    def triple_exponential_func(self, x,a,b,c,d,e,f,g,h,i,j):
        return a*np.exp(-b*(x-c))+d*np.exp(-e*(x-f))+g*np.exp(-h*(x-i))+j

    def ShapeGaussian(self, data, sigma):
        gaussian = []
        for x in self.Time:
            gaussian.append(np.exp(-((x/sigma)**2)/2)*(1/np.sqrt(2*np.pi*sigma*sigma)))
        result = np.convolve(data, gaussian, mode="same")
        return result.tolist()

    def GetPeak(self, data, pol):
        if(pol==1):
            self.Peak = max(data)
        else:
            self.Peak = min(data)
        self.PeakTime = self.Time[data.index(self.Peak)]
        return self.Peak, self.PeakTime

    def GetPeakTime(self, data, pol):
        if(pol==1):
            self.Peak = max(data)
        else:
            self.Peak = min(data)
        self.PeakTime = self.Time[data.index(self.Peak)]
        return self.PeakTime

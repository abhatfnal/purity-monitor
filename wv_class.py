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
    def __init__(self, Files, Pol, Label, ID, VScale = "m", TScale = "u"):
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
        self.Plot = False
        self.ID = ID

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
        self.BaseCounts = self.Samples/10
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
        for x,i in enumerate(range(self.Samples)):
            self.MeanAmp.append(np.average([self.Amp[j][i] for j in range(len(self.Amp))]))
        _,_ = self.GetPeak(self.MeanAmp, self.Pol)
        # PltWfm(time=self.Time, data=self.MeanAmp, label='Signal', xlabel='Time [$\mu$s]', ylabel='Amplitude [mV]',xlim=1,xlim2=1,ylim=-10,ylim2=10)

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
        # f.SetParLimits(0,-1000,1000);
        f.SetParName(0, "Amp");
        f.SetParameter(1,.001);
        # f.SetParLimits(1,-10,10);
        f.SetParName(1, "Tau");
        f.SetParameter(2,0);
        # f.SetParLimits(2,-100,100);
        f.SetParName(2, "XOff");
        f.SetParameter(3,0);
        # f.SetParLimits(3,-100,100);
        f.SetParName(3, "YOff");

        for i in range(repeats):
            hist.Fit("f","RME","");
        fit = self.exponential_func(np.asarray(self.Time), f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3))
        if(self.Plot):
            PPltWfm(self.Time, data, fit, 'Cathode', 'Fit','Time [$\mu$s]', 'Amplitude [mV]',scale=1.2,xlim=self.Time[0],xlim2=self.Time[-1],ylim=min(data)*1.2,ylim2=max(data)*1.2)
        return fit

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

    def FitFullCurve(self, data, start, end, repeats, state=False):
        if(state): print " | Fitting function to total curve..."
        c1 = ROOT.TCanvas( 'c1', 'The Fit Canvas', 200, 10, 700, 500 )
        hist = ROOT.TH1F("waveform", "waveform", self.Samples, self.Time[0],self.Time[-1])
        for i in range(len(data)):
            hist.AddBinContent(i, data[i])
        hist.Draw()


        c1.Update()
        f = ROOT.TF1("f", "([3]/2)*exp(0.5*([1]/[2])^2 - (x-[4])/[2]) * erfc(([1]/[2]-(x-[4])/[1])/sqrt(2))+[0]", start, end)
        f.SetLineColor(ROOT.kRed)
        # f.SetNpx(10000)
    	# f.SetNumberFitPoints(10000);
        f.SetParameter(0, 0.0)
        # f.SetParLimits(0,-1, 1)
        f.SetParName(0, "Baseline")

        f.SetParameter(1,1)
        f.SetParLimits(1,0, 8)
        f.SetParName(1, "Rise Time")

        f.SetParameter(2, 100)
        f.SetParLimits(2, 100, 100000)
        f.SetParName(2, "Fall Time")

        if(np.abs(np.min(data)) > np.abs(np.max(data))):
            f.SetParameter(3, np.min(data))
            f.SetParLimits(3, np.min(data)*1.1, np.min(data)*0.9)
        else:
            f.SetParameter(3, np.max(data))
            f.SetParLimits(3, np.max(data)*0.9, np.max(data)*1.1)
        f.SetParName(3, "Amplitude")

        f.SetParameter(4, 10)
        f.SetParLimits(4, 0, 100000)
        f.SetParName(4, "Peak time")

        for i in range(repeats):
            sys.stdout.write(" | Fit repition... %d" % i)
            sys.stdout.flush()
            sys.stdout.write("\b")
            hist.Fit("f", "REQM", "")
            c1.Update()
        hist.Fit("f", "REM", "")

        f.Draw('SAME')
        c1.Update()
        c1.SaveAs("test.pdf")
        print ' | Extremum:', f.GetMaximum(), f.GetMinimum()
        print ' | Position:', f.GetMaximumX(), f.GetMinimumX()
        print ' | Graph printed...Press any key to continue...'
        raw_input()
        # ROOT.gBenchmark.Show( 'fit1' )
        print " | Reduced chi square...", f.GetChisquare()/f.GetNDF()
        print " | Fit parameters...", f.GetParameter(0), "\t", f.GetParameter(1), "\t", f.GetParameter(2), "\t", f.GetParameter(3), "\t", f.GetParameter(4), "\t", f.GetChisquare()/f.GetNDF()
        # fit = self.FuncExpGausMulti(np.asarray(self.Time), f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3), f.GetParameter(4))
        # if(self.Plot):
        #     PPltWfm(self.Time, data, fit, 'Data', 'Fit','Time [$\mu$s]', 'Amplitude [mV]', scale=1.2, xlim=self.Time[0], xlim2=self.Time[-1], ylim=min(data)*1.2, ylim2=max(data)*1.2, save=True)

        return f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3), f.GetParameter(4)

    def FitFullCurveDouble(self, data, start, end, repeats, state=False):
        if(state): print " | Fitting function to total curve..."
        c2 = ROOT.TCanvas( 'c2', 'The Fit Canvas', 200, 10, 700, 500 )
        hist = ROOT.TH1F("waveform2", "waveform2", self.Samples, self.Time[0],self.Time[-1])
        for i in range(len(data)):
            hist.AddBinContent(i, data[i])

        # p0,p1,p2,p3,p4 = self.FitFullCurve(data, -450, 80, repeats, state=False)
        # p5,p6,p7,p8,p9 = self.FitFullCurve(data, 120, 450, repeats, state=False)


        p0,p1,p2,p3,p4 = 0.0, 2.44765e+00, 2.95234e+01, -4.97247e+01, 1.08153e+01
        p5,p6,p7,p8,p9 = 0.0, 4.49133e+00, 1.33205e+02, 9.92916e+00,  1.21494e+02



        f = ROOT.TF1("f", "[3]/2*exp(0.5*([1]/[2])^2 - (x-[4])/[2]) * erfc(([1]/[2]-(x-[4])/[1])/sqrt(2))+[0] + [8]/2*exp(0.5*([6]/[7])^2 - (x-[9])/[7]) * erfc(([6]/[7]-(x-[9])/[6])/sqrt(2))+[5]", start, end)
        f.SetLineColor(ROOT.kRed)
        f.SetNpx(100000)

        f.FixParameter(0, p0)
        # f.SetParLimits(0,-1, 1)
        f.SetParName(0, "Baseline")

        f.SetParameter(1, p1)
        # f.SetParLimits(1,0, 20)
        f.SetParName(1, "Rise Time")

        f.SetParameter(2, p2)
        # f.SetParLimits(2, 5, 40)
        f.SetParName(2, "Fall Time")

        f.SetParameter(3, p3)
        # f.SetParLimits(3, -5, 0)
        f.SetParName(3, "Amplitude")

        f.SetParameter(4, p4)
        # f.SetParLimits(4, 0, 90)
        f.SetParName(4, "Peak time")

        f.FixParameter(5, p5)
        # f.SetParLimits(0,-1, 1)
        f.SetParName(5, "Baseline 2")

        f.FixParameter(6, p6)
        # f.SetParLimits(1,0, 20)
        f.SetParName(6, "Rise Time 2")

        f.SetParameter(7, p7)
        # f.SetParLimits(2, 5, 40)
        f.SetParName(7, "Fall Time 2")

        f.SetParameter(8, p8)
        # f.SetParLimits(3, -5, 0)
        f.SetParName(8, "Amplitude 2")

        f.SetParameter(9, p9)
        # f.SetParLimits(4, 0, 90)
        f.SetParName(9, "Peak time 2")

        for i in range(repeats):
            print i
            hist.Fit("f", "REQ", "")
        hist.Fit("f", "RE", "")
        hist.Draw()
        f.Draw('SAME')
        c2.Update()
        raw_input()
        print f.GetChisquare(), f.GetNDF(), f.GetChisquare()/f.GetNDF()
        # fit = self.FuncExpGausMulti(np.asarray(self.Time), f.GetParameter(0), f.GetParameter(1), f.GetParameter(2), f.GetParameter(3), f.GetParameter(4))
        if(self.Plot):
            PPltWfm(self.Time, data, fit, 'Data', 'Fit','Time [$\mu$s]', 'Amplitude [mV]', scale=1.2, xlim=self.Time[0], xlim2=self.Time[-1], ylim=min(data)*1.2, ylim2=max(data)*1.2, save=True)

    def FuncExpGausMulti(self, data, p0, p1, p2, p3, p4):
    	#This produces the function 3.5
    	# // p0: baseline, B
    	# // p1: gaussian sig
    	# // p2: exponential decay constant
    	# // p3: number of pulses
    	# // p[4+2*i]: pulse[i] amplitude

    	val = np.full((data.size),p0)
    	time = np.asarray([x-p4 for x in data])
    	#This is the term time = (t - mu)/tau
    	val_tot = val + p3/2.* np.exp((p1*p1/p2/p2/2.-time)/p2)* spl.erfc((p1/p2-time/p1)/np.sqrt(2))+p0
                #This is the function is defined in eqn 3.5
    	return val_tot

    def DoubleFuncExpGausMulti(self, data, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9):
    	#This produces the function 3.5
    	# // p0: baseline, B
    	# // p1: gaussian sig
    	# // p2: exponential decay constant
    	# // p3: number of pulses
    	# // p[4+2*i]: pulse[i] amplitude

    	val = np.full((data.size),p0)
    	time = np.asarray([x-p4 for x in data])
    	#This is the term time = (t - mu)/tau
    	val_tot = val + p3/2.* np.exp((p1*p1/p2/p2/2.-time)/p2)* spl.erfc((p1/p2-time/p1)/np.sqrt(2))+p0
                #This is the function is defined in eqn 3.5
    	return val_tot

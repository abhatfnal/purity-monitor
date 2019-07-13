import numpy as np


def GetFitAmp(self, data, num):
    ScaleUp = 2.0
    ScaleLow = 0.1
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
    import ROOT
    # ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(111)


    if(state): print " | Fitting function to total curve..."
    hist = ROOT.TH1F("waveform", "waveform", self.Samples, self.Time[0],self.Time[-1])
    for i in range(len(data)):
        hist.AddBinContent(i, data[i])
    f = ROOT.TF1("f", "([3]/2)*exp(0.5*([1]/[2])^2 -(x-[4])/[2])*erfc(([1]/[2]-(x-[4])/[1])/sqrt(2))+[0]", start, end)
    f.SetLineColor(ROOT.kRed)
    f.SetNpx(10000)
    f.SetNumberFitPoints(10000)

    ParNames = ['Baseline', 'Rise Time', 'Fall Time', 'Amplitude', 'Peak Time']
    ParValues = [0.0, 1.0, 100.0, GetFitAmp(self, data,0), 10.0]
    ParLow = [-1, 0, 10.0, GetFitAmp(self, data,1), 0]
    ParHigh = [1, 25, 100000.0, GetFitAmp(self, data,2), 1000]
    
    print ParValues 
    print ParLow 
    print ParHigh

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
    import ROOT
    # ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(111)

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



def FitSine(self, data, start, end, repeats, state=False):
    import ROOT
    # ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(111)

    if(state): print " | Fitting function to total curve..."
    hist = ROOT.TH1F("waveform", "waveform", self.Samples, self.Time[0],self.Time[-1])
    for i in range(len(data)):
        hist.AddBinContent(i, data[i])
    f = ROOT.TF1("f", "[0]*sin(x*[1])+[2]", start, end)
    f.SetLineColor(ROOT.kRed)
    f.SetNpx(10000)
    f.SetNumberFitPoints(10000);

    ParNames = ['Amplitude', 'Frequency', 'Baseline']
    ParValues = [1.0, 1.0, 0]
    ParLow = [0, 0.0, -100.0]
    ParHigh = [10000, 150000, 100.0]

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

def FitExponential(self, data, start, end, repeats, state=False):
    import ROOT
    # ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetOptTitle(111)

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
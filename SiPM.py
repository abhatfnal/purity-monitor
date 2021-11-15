import numpy as np
import glob
import h5py

from scipy.signal import butter, lfilter, freqz, filtfilt
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from WaveformAnalysis import Dataset

class SiPM(Dataset.Dataset):
    def __init__(self, Path, Selection='*'):
        self.Path = Path
        self.Selection = Selection
        self.Files = glob.glob(self.Path+self.Selection)
        self.fit_parameters = []
        self.fit_covariance = []
        self.peak_pos = []
        self.peak_height = []
        self.max = []

    def get_filtered_waveform(self, time, amp, lowfreq=1, highfreq=10000, order=6, type='low'):
        fs = 1/(np.mean(np.diff(time))/self.Ch[0].TScale)
        return self.butter_filter(amp, fs, lowfreq, highfreq, order=order, type=type)

    def get_averaged_waveform(self, time, amp, avg=4):
        avg_time = time[:-avg+1]
        avg_amp = np.convolve(amp,np.ones(avg),'valid')/avg
        return avg_time, avg_amp

    def run_fit(self, time, amp):
        max_val = np.max(Amp[cut])
        max_pos_cut = np.where(Amp[cut] == np.max(Amp[cut]))[0][0]
        max_pos = D0.Ch[1].Time[cut][max_pos_cut]

    def butter_filter(self, data, fs, lowfreq, highfreq, order=6, type='low'):
        nyq = 0.5 * fs
        b, a = butter(order, [lowfreq/nyq, highfreq/nyq], btype=type, analog=False)
        y = lfilter(b, a, data)
        return y
    
    def fit_peaks(self, time, data):
        peaks,pdict = find_peaks(data, height=50, width=20, distance=50)
        self.peak_pos.append(peaks)
        self.peak_height.append(pdict['peak_heights'])
        
        for ii,p in enumerate(peaks): 
            pp = time[p]
            ph = pdict['peak_heights'][ii]
            fit = np.where((time>time[p]-50) & (time<time[p]+50))
            try:
                popt, pcov = curve_fit(self.func, time[fit], data[fit], p0=[0, ph, 1, 40, pp], maxfev=10000)
                self.fit_parameters.append(popt)
                self.fit_covariance.append(pcov)
                self.max.append( np.max(self.func(time[fit], *popt)) )
            except:
                self.fit_parameters.append([0,0,0,0,0])
                self.fit_covariance.append([])
    
    def func(self,x,base,V0,sigma,tau,mu):
        return base + V0/2 * np.exp(0.5 * (sigma/tau)**2 - (x-mu)/tau) * erfc(1/np.sqrt(2) * (sigma/tau - (x-mu)/sigma))
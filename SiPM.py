import numpy as np
import glob
import h5py

from scipy.signal import butter, lfilter, freqz, filtfilt
from scipy.special import erfc

from WaveformAnalysis import Dataset


class SiPM(Dataset.Dataset):
    def __init__(self, Path, Selection='*'):
        self.Path = Path
        self.Selection = Selection
        self.Files = glob.glob(self.Path+self.Selection)

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
    
    def func(x,V0,sigma,tau,mu):
        return V0 * np.exp(0.5 * (sigma/tau)**2 - (x-mu)/tau) * erfc(1/np.sqrt(2) * (sigma/tau - (x-mu)/sigma))
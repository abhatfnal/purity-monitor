import numpy as np
import time
import datetime
import sys
import os
import glob
import struct
import matplotlib.pyplot as plt
from optparse import OptionParser
from scipy.signal import blackman
import scipy
from scipy import signal
import plots
from scipy import integrate

from scipy.fftpack import rfft, irfft, fftfreq, fft
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter, sosfilt

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
(options, args) = parser.parse_args()


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sosfilt(sos, data)
    return y

def exponential_func(x, a, b, c, d):
    return a*np.exp(-b*(x-c))+d

def double_exponential_func(x, a, b,  d, e, g):
    return a*np.exp(-b*(x))+d*np.exp(-e*(x))+g

def triple_exponential_func(x,a,b,c,d,e,f,g,h,i,j):
    return a*np.exp(-b*(x-c))+d*np.exp(-e*(x-f))+g*np.exp(-h*(x-i))+j

def GetFFTWaveform(data, lowcut, highcut, sampling, points, order):
    T = points/sampling
    time = np.linspace(-T, T, points, endpoint=False)
    return butter_bandpass_filter(data, lowcut, highcut, sampling, order=order).tolist()

os.system("ls "+options.filepath+"*.txt > list.dat")
num_lines = sum(1 for line in open("list.dat"))
print "Number of files: ", num_lines
file = open("list.dat","r")

time = []
ch1 = []
ch2 = []

for i,line in enumerate(file):
    wfm = open(line[:-1])
    num_points = sum(1 for line in open(line[:-1]))
    for j,data in enumerate(wfm):
        columns = data.split()
        if(i==0):
            time.append(float(columns[0]))
            ch1.append(float(columns[1]))
            ch2.append(float(columns[2]))
        else:
            ch1[j]=ch1[j]+float(columns[1])
            ch2[j]=ch2[j]+float(columns[2])
    wfm.close()
file.close()

time = [x*1000000 for x in time]
baseline_pts = time.index(-100) #baseline until specified time in ms
ch1_b = np.average(ch1[:baseline_pts])/num_lines
ch2_b = np.average(ch2[:baseline_pts])/num_lines

ch1 = [(x/num_lines - ch1_b)*1000 for x in ch1]
ch2 = [(x/num_lines - ch2_b)*1000 for x in ch2]
print "Ch1 Baseline and RMS:    ", ch1_b, np.std(ch1[:baseline_pts])
print "Ch1 Baseline and RMS:    ", ch1_b, np.std(ch1[:baseline_pts])
# ch2 = [0 if (i>time.index(0) and i<(next(ch2.index(i) for i in ch2[time.index(40):] if (i>0 and i<1)))) else x for i,x in enumerate(ch2)]


# plots.plot_single_waveform(time, ch1, 'Time [ms]', 'Amplitude [mV]', 'Cathode', 400, 200)
plots.plot_double_waveform(time, ch1, ch2, 'Cathode', 'Anode','Time [ms]', 'Amplitude [mV]',  200, 200)




N = num_points
T = 1.0 / 5000000.0
x = np.linspace(0.0, N*T, N)
yf = fft(ch1)
yf2 = fft(ch2)
w = blackman(N)
ywf = fft(ch1*w)
ywf2 = fft(ch2*w)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)


fs = 5000000.0
lowcut = 0
highcut = 100E3
T = 0.002
nsamples = T * fs
t = np.linspace(-T, T, nsamples, endpoint=False)

y = GetFFTWaveform(ch1,lowcut=0,highcut=200E3,sampling=5000000.0,points=num_points,order=6)
y2 = GetFFTWaveform(ch2,lowcut=0,highcut=200E3,sampling=5000000.0,points=num_points,order=6)
yf3 = fft(y)
yf4= fft(y2)

print "Ch1 Baseline and RMS FFT:    ", np.average(y[:baseline_pts]), np.std(y[:baseline_pts])
print "Ch2 Baseline and RMS FFT:    ", np.average(y2[:baseline_pts]), np.std(y2[:baseline_pts])

start = -100
end = 15

start1 = 100
end1 = 900

start2 = time[y.index(min(y))]
end2 = 900

start3 = time[y2.index(max(y2))]
end3 = 750



# popt, pcov = curve_fit(exponential_func, np.asarray(time[time.index(start):time.index(end)]), y[time.index(start):time.index(end)], p0=(-1, 1e-6, 0, 1))

# popt1, pcov1 = curve_fit(exponential_func, np.asarray(time[time.index(start1):time.index(end1)]), y[time.index(start1):time.index(end1)], p0=(-.1, 1e-4, 1, 1))
#
# popt2, pcov2 = curve_fit(double_exponential_func, np.asarray(time[time.index(start2):time.index(end2)]), y[time.index(start2):time.index(end2)], p0=(-1.0, 0.1, -1.0, 0.1, 1.0))
#
# popt3, pcov3 = curve_fit(exponential_func, np.asarray(time[time.index(start3):time.index(end3)]), y2[time.index(start3):time.index(end3)], p0=(1, 1e-4, 0, 1))


# popt, pcov = curve_fit(triple_exponenial_func, np.asarray(time[time.index(start):time.index(end)]), y[time.index(start):time.index(end)], p0=(1, 1e-6, 1, 1, 1e-6,1, 1e-6, 1,1,1))
#
#
#



# yy = exponential_func(np.asarray(time), *popt)
# yy1 = exponential_func(np.asarray(time), *popt1)
#
# yy2 = double_exponential_func(np.asarray(time), *popt2)
# yy3 = exponential_func(np.asarray(time), *popt3)


# print "Parameters:  ", popt
# print "Parameters2:  ", popt2
# print "Cathode Value Double at 0:  ", double_exponential_func(time[y.index(min(y))+1],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4])
# idx = (np.abs(y-value)).argmin()
# print "Cathode Value at 0:\t", exponential_func(time[y.index(min(y))],popt1[0],popt1[1],popt1[2],popt1[3]), '\t', min(y), '\t', time[y.index(min(y))]
# print "Anode Value at 0:\t", exponential_func(time[y2.index(max(y2))],popt3[0],popt3[1],popt3[2],popt3[3]), max(y2), '\t', time[y2.index(max(y2))]



# sum = [x+z for x,z in zip(y,y2)]
print len(t), len(ch1)

fig = plt.figure(figsize=(12,7))
ax = fig.gca()
ax.grid()
plt.xlabel("Time [us]")
plt.ylabel("Voltage [mV]")
plt.xlim(-200,200)
plt.ylim(-100,100)
plt.plot(time, y)
plt.plot(time, y2)
# plt.plot(time, yy)
# plt.plot(time, yy3)
# plt.plot(time, yy1)
# plt.plot(time, yy2)
plt.legend(['Cathode', 'Anode','Falling Anode','Falling Cathode', 'Falling Cathode Double'])
plt.show()




fig = plt.figure(figsize=(12,7))

plt.subplot(421)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, ch1, '-b')

plt.subplot(423)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, ch2, '-r')

plt.subplot(422)
plt.grid(True)
plt.title("Frequency band pass filter between 0 and 100 kHz")
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf[1:N//2]), '-b')

plt.subplot(424)
plt.grid(True)
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf2[1:N//2]), '-r')

plt.subplot(425)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, y, '-b')

plt.subplot(427)
plt.grid(True)
plt.xlabel("Time [ms]")
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3,y2, '-r')

plt.subplot(426)
plt.grid(True)
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf3[1:N//2]), '-b')

plt.subplot(428)
plt.grid(True)
plt.xlabel("Frequency [kHz]")
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf4[1:N//2]), '-r')

plt.show()

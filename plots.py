import matplotlib.pyplot as plt
import numpy as np


def plot_single_waveform(time, data, label, xlabel, ylabel, xlim, xlim2, ylim, ylim2):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim,xlim2)
    plt.ylim(ylim,ylim2)
    plt.plot(time, data)
    plt.legend([label])
    plt.show()

def PPltWfm(time, data, data2, label, label2, xlabel, ylabel, scale=1.2, xlim=1, xlim2=1, ylim=1, ylim2=1):
    if(xlim==1 and xlim2==1 and ylim==1 and ylim2==1):
        xlim = time[0]
        xlim2 = time[-1]
        ylim = min(min(data),min(data2))*scale
        ylim2 = max(max(data),max(data2))*scale
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim,xlim2)
    plt.ylim(ylim,ylim2)
    plt.plot(time, data)
    plt.plot(time, data2)
    plt.legend([label, label2])
    plt.show()

def plot_single_fft(time, data, data2, label, label2, xlabel, ylabel, xlim, xlim2, ylim, ylim2):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim,xlim2)
    plt.ylim(ylim,ylim2)
    plt.loglog(time, data)
    plt.loglog(time, data2)
    plt.legend(['Cathode', 'Anode','Sum'])
    plt.show()

def plot_double_fft(time, freq, data, data2, fft, fft2, label, label2, xlabel, ylabel, xlim, ylim):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(-xlim,xlim)
    plt.ylim(-ylim,ylim)
    plt.plot(time, data)
    plt.plot(time, data2)
    plt.legend(['Cathode', 'Anode','Sum'])

    plt.subplot(421)
    plt.ylabel("Amplitude [mV]")
    plt.plot(t*1E3, ch1, '-b')

    plt.subplot(423)
    plt.ylabel("Amplitude [mV]")
    plt.plot(t*1E3, ch2, '-r')

    plt.subplot(422)
    plt.title("Frequency band pass filter between 0 and 100 kHz")
    plt.ylabel("Amplitude")
    plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf[1:N//2]), '-b')

    plt.subplot(424)
    plt.ylabel("Amplitude")
    plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf2[1:N//2]), '-r')

    plt.subplot(425)
    plt.ylabel("Amplitude [mV]")
    plt.plot(t*1E3, y, '-b')

    plt.subplot(427)
    plt.xlabel("Time [ms]")
    plt.ylabel("Amplitude [mV]")
    plt.plot(t*1E3,y2, '-r')

    plt.subplot(426)
    plt.ylabel("Amplitude")
    plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf3[1:N//2]), '-b')

    plt.subplot(428)
    plt.grid(True)
    plt.xlabel("Frequency [kHz]")
    plt.ylabel("Amplitude")
    plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf4[1:N//2]), '-r')

    plt.show()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


def PltScatter(xvalue, yvalue, label, xlabel, ylabel, scale=1.2, xlim=1, xlim2=1, ylim=1, ylim2=1, save=False):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.scatter(xvalue, yvalue)
    plt.legend([label])

def PltScatterD(xvalue, yvalue, yvalue2, label, xlabel, ylabel, scale=1.2, xlim=1, xlim2=1, ylim=1, ylim2=1, save=False):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    plt.xlabel('Time [$\mu$s]')
    plt.ylabel('Amplitude [mV]')
    plt.scatter(xvalue, yvalue)
    plt.scatter(xvalue, yvalue2, color='red')
    plt.legend([label])

def PltWfm(time, data, legend, scale=1.5, xlim=1, xlim2=1, ylim=1, ylim2=1, save=False):
    if(xlim==1 and xlim2==1):
        xlim = time[0]
        xlim2 = time[-1]
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    plt.xlabel('Time [$\mu$s]')
    plt.ylabel('Amplitude [mV]')
    plt.xlim(xlim,xlim2)
    for signal in data:
        plt.plot(time, signal)
    plt.legend(legend)

def PPltWfm(time, data, data2, label, label2, xlabel, ylabel, scale=1.2, xlim=1, xlim2=1, ylim=1, ylim2=1, save=False):
    if(xlim==1 and xlim2==1 and ylim==1 and ylim2==1):
        xlim = time[0]
        xlim2 = time[-1]
        ylim = min(min(data),min(data2))*scale
        ylim2 = max(max(data),max(data2))*scale
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim,xlim2)
    plt.ylim(ylim,ylim2)
    plt.plot(time, data, color='b', linewidth='1')
    plt.plot(time, data2, color='r', linewidth='1')
    plt.legend([label, label2])
    # plt.show()
    # plt.savefig("double.pdf")

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
    # plt.show()

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

    # plt.show()

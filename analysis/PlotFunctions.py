import matplotlib, os, datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SavePath = '/home/fas/david_moore/aj487/purity_monitor/plots/analysis/'
Date = datetime.datetime.now().strftime("%Y%m%d")

# plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex':True,'font.size':12,'font.family':'serif'}
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

def PltTime(Time,Data,Legend,Label,XRange=0,YRange=0,XTicks=0,YTicks=0,SaveName='',Save=False):
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=XTicks))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(YTicks))

    ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.7)
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    plt.xlabel('Time [hh:mm:ss]', fontsize=14)
    plt.ylabel(Label, fontsize=14)

    plt.xlim(np.min(Time),np.max(Time))
    formatter = matplotlib.dates.DateFormatter('%H:%M:%S')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.gcf().autofmt_xdate()
    if YRange != 0: 
        plt.ylim(YRange[0], YRange[1])
    else:
        maxY = np.max([np.max(x) for x in Data])
        plt.ylim(0, maxY*1.2)
    
    for ii, yy in enumerate(Data):
        plt.plot(Time, yy, 'o', label=Legend[ii], color=colors[ii], marker='o', mew=0.01, markersize=4)
    plt.legend(loc='lower left',bbox_to_anchor=(0.005,0.92), ncol=3, borderaxespad=0, fontsize=14)
    fig.tight_layout()
    if Save:
        plt.savefig(SavePath+Date+'/'+SaveName+'.pdf',bbox_inches='tight')

def PltScatter(xvalue, yvalue, legend, xlabel, ylabel, scale=1.2, xlim=1, xlim2=1, ylim=1, ylim2=1, save=False):
    fig = plt.figure(figsize=(6,5))
    ax = fig.gca()
    ax.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    color = ['red', 'blue']
    for signal in yvalue:
        plt.scatter(xvalue, signal)
    plt.legend(legend)
    fig.tight_layout()

def PltWfm(Time,Data,Legend,Label=['Time [$\mu$s]', 'Amplitude [mV]'],XRange=0,YRange=0,XTicks=0,YTicks=0,Color='',SaveName='',Save=False):
    fig = plt.figure(figsize=(10,6))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='-')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(XTicks))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(YTicks))
    plt.xlabel(Label[0])
    plt.ylabel(Label[1])
    if XRange != 0:
        plt.xlim(XRange[0],XRange[1])
    else: 
        plt.xlim(Time[0],Time[-1])
    if YRange != 0:
        plt.ylim(YRange[0],YRange[1])
    for ii,signal in enumerate(Data):
        if Color == '':
            plt.plot(Time, signal, label=Legend[ii], color=colors[ii],linewidth=2.0)
        if Color == 'k': 
            plt.plot(Time, signal, label='', color='black',linewidth=1.0)
    ax.legend(loc='upper right',fontsize=16)
    fig.tight_layout()
    plt.margins(0,0)
    if Save:
        plt.savefig(SavePath+Date+'/'+Save+'.pdf',bbox_inches='tight')

def PltAllWfm(Time,Data,Legend,Label=['Time [$\mu$s]', 'Amplitude [a.u.]'],XRange=0,YRange=0,XTicks=0,YTicks=0,Color='',Save=''):
    fig = plt.figure(figsize=(16,9))
    ax = fig.gca()
    plt.xlabel(Label[0])
    plt.ylabel(Label[1])
    plt.xlim(Time[0],Time[len(Time)//2])
    offset = 10 
    plt.ylim(1E-1,len(Data)*offset)
    plt.yscale('log')
    for ii,signal in enumerate(Data):
        signal = signal +ii*offset
        signal = np.roll(signal, len(Time)//2 + 100*ii)
        if Color == '':
            plt.plot(Time, signal, color=colors[ii],linewidth=2.0)
        if Color == 'k': 
            plt.plot(Time, signal, color='black',linewidth=.1)
    plt.savefig('wfm_all.png')

def plot_single_fft(time, data, data2):
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid()
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # plt.xlim(xlim,xlim2)
    # plt.ylim(ylim,ylim2)
    plt.loglog(time, data)
    plt.loglog(time, data2)
    plt.legend(['Cathode', 'Anode','Sum'])

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

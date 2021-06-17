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

def PltTime(Time,Data,Legend,Label,XRange=0,YRange=0,XTicks=0,YTicks=0,SaveName='',Title='',Save=False):
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=XTicks))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(YTicks))

    ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.7)
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    plt.xlabel('Time [hh:mm:ss]')
    plt.ylabel(Label)

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
        window = 10
        plt.plot(Time[(window-1):], np.convolve(yy, np.ones(window), 'valid') / window, label=Legend[ii], color=colors[ii], mew=0.01, markersize=4)
        plt.plot(Time, yy, 'o', color=colors[ii], marker='o', mew=0.01, markersize=2, alpha=0.5)
    plt.legend(loc='lower left',bbox_to_anchor=(0.005,0.92), ncol=3, borderaxespad=0)
    plt.title(Title)
    fig.tight_layout()
    if Save:
        plt.savefig(SavePath+Date+'/'+SaveName+'.pdf',bbox_inches='tight')

def PltChargeVsTime(Time,Data,CC,Legend,Label,XRange=0,YRange=0,XTicks=0,YTicks=0,Title=''):
    fig,ax = plt.subplots()

    window = 10
    
    ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.5)
    ax.grid(b=True, which='minor', color='k', linestyle=':', alpha=0.5)
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(YTicks))
    plt.xlabel('Time [hh:mm:ss]')
    plt.ylabel(Label)

    plt.xlim(np.min(Time),np.max(Time))
    formatter = matplotlib.dates.DateFormatter('%H:%M')
    # ax.xaxis.set_major_formatter(formatter)
    # plt.gcf().autofmt_xdate()
    if YRange != 0: 
        plt.ylim(YRange[0], YRange[1])
    else:
        maxY = np.max([np.max(x) for x in Data])
        plt.ylim(0, maxY*1.2)
    
    for ii, yy in enumerate(Data):
        plt.plot(Time, yy, 'o', color=colors[ii], marker='o', mew=0.01, markersize=2, alpha=0.8)
        plt.plot(Time[int((window-1)/2):-int((window-1)/2)-1], np.convolve(yy, np.ones(window), 'valid') / window, label=Legend[ii], color=colors[ii], mew=0.01, markersize=4)
    
    
    plt.title(Title)

    ax2 = ax.twinx()
    ax2.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax2.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=XTicks))
    ax2.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(YTicks))

    ax2.xaxis.set_major_formatter(formatter)
    plt.gcf().autofmt_xdate()
    plt.ylim(0,200)
    plt.ylabel('Charge Collection [\%]')
    ax2.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(20))
    plt.plot(Time, CC, 'o', color=colors[2], marker='o', mew=0.01, markersize=2, alpha=0.8)
    
    run_av = np.convolve(CC, np.ones(window), 'valid') / window
    plt.plot(Time[int((window-1)/2):-int((window-1)/2)-1], run_av, label=Legend[2], color=colors[2], mew=0.01, markersize=4, zorder=0)

    leg = ax.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, borderaxespad=0, framealpha=1.0)
    leg = ax.get_legend()
    leg.remove()
    ax2.add_artist(leg)
    ax2.legend(loc='upper right',bbox_to_anchor=(0.99,0.99), ncol=1, borderaxespad=0, framealpha=1.0).set_zorder(10)

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
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if by_label: 
        ax.legend(loc='upper right')
    fig.tight_layout()
    plt.margins(0,0)
    if Save:
        plt.savefig(SavePath+Date+'/'+SaveName+'.pdf',bbox_inches='tight')

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

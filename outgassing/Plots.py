import numpy as np 
import datetime 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from matplotlib.dates import DateFormatter
from scipy.stats import chisquare
from scipy.optimize import curve_fit

params = {'text.usetex':True,'font.size':12,'font.family':'serif'}
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', '#F27781', '#18298C', '#04BF8A', '#F2CF1D', '#F29F05', '#7155D9', '#8D07F6', '#9E91F2', '#F29B9B', '#F25764', '#6FB7BF', '#B6ECF2', '#5D1314', '#B3640F']

def func(x, a, b, c, d):
    return a * np.exp(-(x-c)/b) +d

def FitExponential(Time, Data, Range):
    Cut = (Time > Range[0]) & (Time < Range[1])
    X = Time[Cut]
    Y = Data[Cut]
    p, pcov = curve_fit(func, X, Y, p0=(1E-5,40, 5, 1E-5), maxfev=1000000)
    return [X, func(X, *p), r'$%.1e\cdot e^{-\frac{t-%.3f}{%.3f}} + %.3e$' % (p[0],p[2],p[1],p[3]), p, pcov]

def PlotBestFitOverTime(Time, Data, XRange=0, YRange=0, Fit=None):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.xlabel('Time [hours]', fontsize=14)
    plt.ylabel(r'Outgassing Rate [mbar$\,\cdot\,$Liter/s]', fontsize=14)
    if XRange == 0:    
        XMin = np.min([np.min(x) for x in Time])
        XMax = np.max([np.max(x) for x in Time])
        plt.xlim(XMin, XMax)
    else: 
        print('test')
        plt.xlim(XRange[0], XRange[1]) 
    if YRange == 0:    
        plt.ylim(1E-10, 1E-5)
    else: 
        plt.ylim(YRange[0], YRange[1]) 

    for jj, (x,y) in enumerate(zip(Time,Data)):
        for ii, key in enumerate(y.keys()): 
            if key not in ['Oxygen', 'Nitrogen','Carbondioxide']: 
                pass  
            else:
                if jj==1:
                    plt.plot(x, y[key], label='Background', color = colors[jj], linewidth=1.0)
                else:
                    plt.plot(x, y[key], label=key, color = colors[ii], linewidth=1.0)
    if Fit is not None:
        plt.plot(Fit[0], Fit[1], label=Fit[2], color = 'red', ls='-')
    plt.legend(loc='upper right', prop={'size': 14}, ncol=1)
    fig.tight_layout()

def PlotRGASpectrum(Time, Mass, Pressure, Dict, ArrowPos): 
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.xlim(0,50)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Mass [u]', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    plt.yscale('log')
    data = Pressure[::400]
    for ii, x in enumerate(data):
        plt.plot(Mass, x, linewidth=0.8, color='%.3f' % float(ii/len(data)))
    for ii, x in enumerate(Dict.keys()):
        ax.annotate(Dict[x]['Formula'], xy=(Dict[x]['Mass'], ArrowPos[x]*1.2), xytext=(Dict[x]['Mass'], ArrowPos[x]*5), 
                    arrowprops=dict(facecolor='black',width=0.6,headwidth=4), fontsize=12, ha='center', 
                    bbox=dict(boxstyle="round", fc="white", ec="grey"))
    # plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)
    fig.tight_layout()

def PlotOverTimeRelative(Time, Data, XRange=0, YRange=0):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.xlabel('Time [hours]', fontsize=14)
    plt.ylabel('Partial Pressure [mbar]', fontsize=14)
    if XRange == 0:    
        XMin = np.min([np.min(x) for x in Time])
        XMax = np.max([np.max(x) for x in Time])
        plt.xlim(XMin, XMax)
    else: 
        plt.xlim(XRange[0], XRange[1]) 
    if YRange == 0:    
        plt.ylim(1E-10, 1E-5)
    else: 
        plt.ylim(YRange[0], YRange[1]) 

    for jj, (x,y) in enumerate(zip(Time,Data)):
        for ii, key in enumerate(y.keys()): 
            if key not in ['Oxygen', 'Nitrogen','Carbondioxide']: 
                pass  
            else:
                if jj==1:
                    plt.plot(x, y[key], label='Background', color = colors[jj], linewidth=1.0)
                else:
                    plt.plot(x, y[key], label=key, color = colors[ii], linewidth=1.0)
    plt.legend(loc='upper right', prop={'size': 14}, ncol=1)
    fig.tight_layout()


def PlotOverTimeAbsolut(Time,Data, XRange=0, YRange=0): 
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    # ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gcf().autofmt_xdate()
    formatter = DateFormatter('%H:%M:%S')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

    plt.xlabel('Time [HH:MM:SS]', fontsize=14)
    plt.ylabel('Partial Pressure [mbar]', fontsize=14)
    # print(Time)
    print(np.shape(Time))
    print(np.shape(Data))
    # quit()
    if XRange == 0:    
        XMin = np.min([np.min(x) for x in Time])
        XMax = np.max([np.max(x) for x in Time])
        plt.xlim(XMin, XMax)
    else: 
        plt.xlim(XRange[0], XRange[1]) 
    if YRange == 0:    
        plt.ylim(1E-10, 1E-5)
    else: 
        plt.ylim(YRange[0], YRange[1]) 

    for jj, (x,y) in enumerate(zip(Time,Data)):
        for ii, key in enumerate(y.keys()): 
            if key not in ['Oxygen', 'Nitrogen','Carbondioxide', 'Water']: 
                pass  
            else:
                if jj==1:
                    plt.plot(x, y[key], label='Background', color = colors[jj], linewidth=1.0)
                else:
                    plt.plot(x, y[key], label=key, color = colors[ii], linewidth=1.5)
    plt.legend(loc='upper right', prop={'size': 14}, ncol=1)
    fig.tight_layout()


    # fig = plt.figure(figsize=(12.4,7))
    # ax = fig.gca()
    # ax.grid(b=True, which='major', color='k', linestyle='--')
    # ax.grid(b=True, which='minor', color='grey', linestyle=':')
    # ax.set_yscale('log')
    
    # plt.xticks(fontsize = 14)
    # plt.yticks(fontsize = 14)
    # plt.xlabel('Time', fontsize=14)
    # plt.ylabel('Pressure [mbar]', fontsize=14)
    # timet = [datetime.datetime.strptime(x, "%Y%m%d%H%M%S") for x in Time]
    # plt.xlim(timet[0], timet[-1])
    # plt.ylim(1E-10, 1E-4)
    # for key in Data.keys(): 
    #     plt.plot(timet, Data[key], label=key)
    # plt.legend(loc='upper left', prop={'size': 14}, ncol=4)
    # fig.tight_layout()

def PlotScatter(X, Y, Labels=[], XRange=[], YRange=[], Legend=[]):
    fig = plt.figure(figsize=(10,8))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel(Labels[0], fontsize=14)
    plt.ylabel(Labels[1], fontsize=14)
    plt.xlim(XRange[0], XRange[1])
    # plt.ylim(YRange[0], YRange[1])
    for ii, (x,y) in enumerate(zip(X,Y)):
        plt.scatter(x,y, label=Legend[ii], color=colors[ii])

    plt.legend(loc='upper left', prop={'size': 14}, ncol=4)
    fig.tight_layout()
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

def func(x, a, b, c):
    return a * np.exp(-x/b) + c

def PlotBestFitOverTime(Time, Data):
    End = 5000
    fig = plt.figure(figsize=(12.4,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.xlabel('Time [hours]', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    timet = [datetime.datetime.strptime(x, "%Y%m%d%H%M%S") for x in Time]
    timet = np.array([datetime.timedelta.total_seconds(x-timet[0]) for x in timet])
    timet = timet/3600.0

    ttime = timet[:End]
    plt.xlim(ttime[0], ttime[-1])
    plt.ylim(1E-8, 1E-6)
    
    for ii, key in enumerate(Data.keys()): 
        if key is not 'Oxygen': 
            pass  
        else:
            data = Data[key][:End]
            plt.plot(ttime, data, label=key, color = colors[ii], linewidth=1.0)
            p, pcov = curve_fit(func, ttime, data, p0=(1E-8,1,1E-8), maxfev=1000000)
            xsquared = chisquare(data,func(ttime,*p), axis = None)
            # print(str(label3[ii]) + str(param) + " chi squared value = " + str(xsquared) + "\n")
            plt.plot(ttime, func(ttime, *p), label=r'$%.1e\cdot\exp(-t/%.2f)+%.1e$' % (p[0],p[1],p[2]), color = colors[ii], linestyle='--')
    plt.legend(loc='upper right', prop={'size': 14}, ncol=1)


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


def PlotOverTime(Time, Data):
    fig = plt.figure(figsize=(12.4,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.set_yscale('log')
    plt.gcf().autofmt_xdate()
    formatter = DateFormatter('%b %d %H:%M')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Pressure [mbar]', fontsize=14)
    timet = [datetime.datetime.strptime(x, "%Y%m%d%H%M%S") for x in Time]
    plt.xlim(timet[0], timet[-1])
    plt.ylim(1E-10, 1E-4)
    for key in Data.keys(): 
        plt.plot(timet, Data[key], label=key)
    plt.legend(loc='upper left', prop={'size': 14}, ncol=4)
    fig.tight_layout()
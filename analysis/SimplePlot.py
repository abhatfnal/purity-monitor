import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
import pandas as pd 
from scipy.optimize import curve_fit
import numpy as np

params = {'text.usetex':True,'font.size':14,'font.family':'serif'}
plt.rcParams.update(params) 

plt.rcParams.update({'errorbar.capsize': 2})

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

def Plot(X,Y,Labels,Legend,Ticks,XLim,YLim): 
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlim(XLim[0],XLim[1])
    plt.ylim(YLim[0],YLim[1])
    plt.xlabel(Labels[0], fontsize=16)
    plt.ylabel(Labels[1], fontsize=16)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(Ticks[0]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Ticks[1]))
    plt.grid()
    for ii,(x,y) in enumerate(zip(X,Y)):
        plt.scatter(x, y, label=Legend[ii], color=colors[ii], marker='d')
    ax.legend(loc='upper right')

def PlotError(X,Y,Err,Labels,Legend,Ticks,XLim,YLim): 
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlim(XLim[0],XLim[1])
    plt.ylim(YLim[0],YLim[1])
    plt.xlabel(Labels[0], fontsize=16)
    plt.ylabel(Labels[1], fontsize=16)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(Ticks[0]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Ticks[1]))
    plt.grid()
    for ii,(x,y) in enumerate(zip(X,Y)):
        plt.errorbar(x, y, yerr=Err[ii], label=Legend[ii], color=colors[ii], marker='o', ms=4)
    ax.legend(loc='upper right')

def PlotFit(X,Y,FitX,FitY,Labels,Legend,Ticks,XLim,YLim): 
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    plt.xlim(XLim[0],XLim[1])
    plt.ylim(YLim[0],YLim[1])
    plt.xlabel(Labels[0], fontsize=16)
    plt.ylabel(Labels[1], fontsize=16)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(Ticks[0]))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Ticks[1]))
    plt.grid()
    for ii,(x,y) in enumerate(zip(X,Y)):
        plt.scatter(x, y, label=Legend[ii], color=colors[ii], marker='d')
    plt.plot(FitX, FitY, label='Fit', color='red')
    ax.legend(loc='upper right')

def ExpNeg(x,a,b):
    return np.exp(-a*x)*b


data = pd.read_csv('$PROJ/purity_monitor/data/20190830.txt', sep=",")
data.columns = ["Voltage", "Drift Time", "Collection", "Lifetime", 'Drift Time Error', 'Collection Error', 'Lifetime Error Low', 'Lifetime Error High']

exo = pd.read_csv('../data/exo200_driftvelocity.txt', sep=",")
exo.columns = ["Electric Field", "Drift Velocity"]

DriftLength = 80 # mm
Resistance = 80 # MOhm
DriftTime = data['Drift Time'] # us
DriftField = data['Voltage'] # total V across DriftLength 
DriftField = (DriftField*Resistance/(200.0+Resistance))/(DriftLength/10.0)
DriftVelocity = DriftLength/DriftTime

data['Drift Time Error'] = [3] * len(data['Drift Time Error'])


Plot(X=[DriftField],Y=[data['Collection']*100],Labels=['Drift Field [V/cm]',r'Charge Collection [\%]'],Legend=['YLXPM'],Ticks=[20,20],XLim=[0,140],YLim=[0,120])
# plt.savefig('collection_vs_driftfield.pdf')

Plot(X=[DriftField],Y=[data['Lifetime']],Labels=['Drift Field [V/cm]',r'Lifetime [$\mu$s]'],Legend=['YLXPM'],Ticks=[20,1000],XLim=[0,140],YLim=[0,15000])
# plt.savefig('lifetime_vs_driftfield.pdf')

popt, cov = curve_fit(ExpNeg, DriftTime, data['Collection']*100, p0=[0.1,100])
print(1.0/popt[0])
x = np.linspace(40,140,100)
PlotFit(X=[DriftTime],Y=[data['Collection']*100],FitX=x,FitY=ExpNeg(x, *popt),Legend=['YLXPM'],Labels=['Drift Time [$\mu$s]',r'Charge Collection [\%]'],Ticks=[20,20],XLim=[20,160],YLim=[0,120])
# plt.savefig('collection_vs_drifttime.pdf')

plt.show()

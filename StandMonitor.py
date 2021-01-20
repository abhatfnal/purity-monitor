import numpy as np
import time, datetime, sys, os, glob, struct
import matplotlib
import matplotlib.pyplot as plt
from optparse import OptionParser
import h5py
import datetime
from matplotlib.dates import DateFormatter
from matplotlib.ticker import MultipleLocator

colors = ['#1f78b4', '#e66101', '#d62728', '#33a02c']

class SensorData:
    def __init__(self, Filepath, PlotTime=0):
        self.Filepath = Filepath
        self.File = h5py.File(self.Filepath, 'r')
        self.Key = list(self.File.keys())[0]
        self.RawData = np.array(self.File[self.Key])
        self.Date = os.path.split(self.Filepath)[1][:-3]
        self.RefTime, self.DateTime = self.GetDateFromInput(self.Date)
        self.Data = {}
        self.Index = np.arange(0,16,1)
        self.PlotTime = PlotTime
        self.StartTime = datetime.datetime.now()
        self.Labels = ['Gas System', 'Chamber', 'Stainless-steel Cylinder 1', 'Stainless-steel Cylinder 2', 'LN Dewar 1', 'LN Dewar 2', 'Xenon Pump', 'Flow Meter', 'Back Pump', 'Cold Head', 'Copper Ring', 'Copper Jacket', 'TPC Bottom', 'dummy1', 'dummy1', 'Time']
        
    def GetData(self, Selection=None):
        if Selection is None: 
            for ii,Label in enumerate(self.Labels):
                if Label == 'Time':
                    self.Seconds = [(x - self.RefTime - 14400 - 3600) for x in self.RawData[:,15]]
                    self.Time = [self.DateTime + datetime.timedelta(seconds=x) for x in self.Seconds]
                else:
                    self.Data[Label] = self.RawData[:,self.Index[ii]]
        else: 
            SelectionIndex = np.where(self.Labels == Selection)[0][0]
            self.Data[Selection] = self.RawData[:,self.Index[SelectionIndex]]

    def GetDateFromInput(self, Date):
        dd = list(Date)
        year = int(dd[0]+dd[1]+dd[2]+dd[3])
        month = int(dd[4]+dd[5])
        day = int(dd[6]+dd[7])
        dt = datetime.datetime(year,month,day,00,00,00)
        tmp = datetime.datetime(1904,1,1,0,0)
        at = int((dt - tmp).total_seconds())
        return at, dt

    def PlotData(self, Selection='Temperature', Time=None, Data=None, XYLabel=None, Labels=None, Tags=None, XRange=0, YRange=[1,1], Ticks=[0,5], XTicks=2):
        if Selection == 'Temperature': 
            XYLabels = ['Time [hh:mm]', 'Temperature [C]']
            Tags = self.Labels[9:13]
            Data = [self.Data[x] for x in Tags]
        elif Selection == 'Xenon Pressure': 
            XYLabels = ['Time [hh:mm]', 'Pressure [PSIG]']
            Tags = self.Labels[2:4]
            Data = [self.Data[x] for x in Tags]
        elif Selection == 'System Pressure': 
            XYLabels = ['Time [hh:mm]', 'Pressure [PSIG]']
            Tags = self.Labels[0:2]
            Data = [self.Data[x] for x in Tags]

        fig = plt.figure(figsize=(12,9))
        ax = fig.gca()
        if(YRange[0]!=1 or YRange[1]!=1):
            plt.ylim(YRange[0], YRange[1])
        
        ax.minorticks_on()
        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=XTicks))
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(Ticks[1]))
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))

        ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.8)
        ax.grid(b=True, which='minor', color='grey', linestyle=':')

        plt.xlabel(XYLabels[0], fontsize=14)
        plt.ylabel(XYLabels[1], fontsize=14)

        plt.gcf().autofmt_xdate()
        formatter = DateFormatter('%H:%M')
        plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

        for ii,(X,Tag) in enumerate(zip(Data,Tags)):
            plt.plot(self.Time, X, label=Tag, linewidth=2, color=colors[ii])
        plt.legend(loc='upper left')

        if(XRange != 0):
            xlim1 = startTime - datetime.timedelta(seconds=60*XRange)
            xlim2 = startTime + datetime.timedelta(seconds=60*XRange/4)
            ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=10))
            ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(interval=5))
        else:
            xlim1 = self.DateTime + datetime.timedelta(seconds=3600*0)
            xlim2 = self.DateTime + datetime.timedelta(seconds=3600*24)
        plt.xlim(xlim1, xlim2)

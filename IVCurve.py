import numpy as np
import scipy
import glob
import h5py
import datetime as dt
import matplotlib.pyplot as plt 

class IVCurve:
    def __init__(self, Path, Selection='*', Go=False):
        self.Path = Path
        self.Selection = Selection
        self.Files = glob.glob(self.Path+self.Selection)
        if Go:
            self.get_data()
            self.get_gradient(self.Voltage,self.Current)
            self.format_timestamp()

    def get_data(self):
        self.Voltage = []
        self.Current = []
        self.Timestamp = []
        for filename in self.Files:
            f = h5py.File(filename, 'r')
            for ii,group in enumerate(f.keys()):
                self.Timestamp.append(group)
                self.Voltage.append(np.array(f[group]['Voltage'][:]))
                self.Current.append(np.array(f[group]['Current'][:]))
            f.close()
    
    def get_gradient(self, Voltage, Current):
        self.Gradient = []
        for x,y in zip(Voltage, Current):
            grad = np.gradient(np.log(y))/np.gradient(x)
            self.Gradient.append(grad)

    def get_average(self, Size=-1):
        if Size == -1:
            Size = len(self.Current)
        self.AvgVoltage = []
        self.AvgCurrent = []
        self.AvgTimestamp = []
#         self.TimeInterval=[]
        for x in range(len(self.Current)):
            ii = x*Size
            if ii >= len(self.Current):
                break
            voltage = np.mean(self.Voltage[ii:ii+Size], axis=0)
            current = np.median(self.Current[ii:ii+Size], axis=0)
            timestamp = np.mean([float(x) for x in self.Timestamp[ii:ii+Size]])
#             timeinterval=np.array([[self.Timestamp[ii]],[self.Timestamp[ii+Size]]])
            self.AvgVoltage.append(voltage)
            self.AvgCurrent.append(current)
            self.AvgTimestamp.append(timestamp)
#             self.TimeInterval.append(timeinterval)
    
    def format_timestamp(self, Format='datetime', Ref=dt.datetime(2021,11,9,0,0)):
        if Format == 'datetime':
            self.Datetime = []
            date = self.Path.split('/')[-2]
            year = int(date[:4])
            month = int(date[4:6])
            day = int(date[-2:])
            for x in self.Timestamp:
                tt = int(x.split('.')[0])
                hour = int(tt/3600.0)
                minute = int((tt-hour*3600)/60)
                second = tt - hour*3600 - minute*60 
                self.Datetime.append(dt.datetime(year,month,day,hour,minute,second))

    def plot(self):
        plt.figure()
        plt.xlabel('Bias Voltage [V]')
        plt.ylabel('Current [A]')
        plt.yscale('log')
        plt.xlim(np.min([np.min(x) for x in self.Voltage]), np.max([np.max(x) for x in self.Voltage]))
        plt.ylim(1e-12,1e-6)
        for x,y in zip(self.Voltage, self.Current):
            plt.plot(x,y, color='k')
        plt.show()
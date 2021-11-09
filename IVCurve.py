import numpy as np
import scipy
import glob
import h5py
import datetime as dt

class IVCurve:
    def __init__(self, Path, Selection='*', Go=False):
        self.Path = Path
        self.Selection = Selection
        self.Files = glob.glob(self.Path+self.Selection)
        if Go:
            self.get_data()
            self.get_gradient(self.Voltage,self.Current)

    def get_data(self):
        self.Voltage = []
        self.Current = []
        self.TimeStamp = []
        for filename in self.Files:
            f = h5py.File(filename, 'r')
            for ii,group in enumerate(f.keys()):
                self.TimeStamp.append(group)
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
        for x in range(len(self.Current)):
            ii = x*Size
            if ii >= len(self.Current):
                break
            voltage = np.mean(self.Voltage[ii:ii+Size], axis=0)
            current = np.median(self.Current[ii:ii+Size], axis=0)
            self.AvgVoltage.append(voltage)
            self.AvgCurrent.append(current)
    
    def format_timestamp(self, Format='datetime', Ref=dt.datetime(2021,11,9,0,0)):
        if Format == 'datetime':
            pass
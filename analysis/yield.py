import numpy as np
import matplotlib.pyplot as plt
import time, datetime, sys, os, glob, struct, h5py, argparse

from WaveformClass import *
from HelperClasses import *
from MatchedFilter import MatchedFilter

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="filepath", nargs="*", help="Specify path/file to be read in.")
parser.add_argument("-d", "--dir", type=str, action="store",  dest="dirpath", nargs="*", help="Analyze directory and subdirectories.")
parser.add_argument("--txt", action="store_true", dest="txt", default=False, help="")
parser.add_argument("-s", action="store_true", dest="show", default=False, help="")
arg = parser.parse_args()

def ImportDataFromHDF5(file, channels):
    f = h5py.File(file, 'r')
    print(" | Filename...", file)
    Keys = list(f.keys())
    for ch in channels:
        ch.Time = np.array(f.get('Time')).flatten() * ch.TScale
        ch.Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
        Group = f.get(ch.ChName)
        GroupKeys = Group.keys()
        # print(" | Number of files in ch%d...\t %d/%d" % (ch.ID, ch.Files[-1], np.sum(ch.Files)))
        for key in GroupKeys:
            ch.Amp.append(np.array(Group.get(key)).flatten() * ch.VScale)
            ch.TimeStamp.append(datetime.datetime.strptime((f.attrs['Date']+Group.get(key).attrs["TimeStamp"]).decode('utf-8'), '%Y%m%d%H%M%S'))
    
def Yield(file): 
    f = h5py.File(file, 'r')
    print(" | Filename...", file)
    Keys = list(f.keys())
    Time = np.array(f.get('Time')).flatten() * ch.TScale
    Trigger = np.array(f.get('Trigger')).flatten() * ch.VScale
    Group = f.get(ch.ChName)
    GroupKeys = Group.keys()
    for key in GroupKeys:
        Amp = np.array(Group.get(key)).flatten() * ch.VScale
        TimeStamp = datetime.datetime.strptime((f.attrs['Date']+Group.get(key).attrs["TimeStamp"]).decode('utf-8'), '%Y%m%d%H%M%S')
        yield (Amp, TimeStamp)
    

def ChooseFilesToAnalyze(arg):
    FileNames = []
    Method = ''
    if (arg.dirpath != None):
        for dir in arg.dirpath:
            FileNames.append(glob.glob(dir+"*.h5"))
        Method = 'Single'
    if (arg.filepath != None):
        if ('*' in arg.filepath):
            FileNames.append(glob.glob(dir+arg.filepath))
        else:
            FileNames.append(arg.filepath)
        Method = 'All'
    FileNames = [val for sublist in FileNames for val in sublist]
    return FileNames, Method

if __name__ == '__main__':

    ###### Initialize class object for each channel.
    ch1 = WFM(1)
    ch2 = WFM(2)
    ch1.Pol = +1
    ch2.Pol = -1
    channels = [ch1, ch2]

    ###### Import data from HDF5 file.
    FileNames, Method = ChooseFilesToAnalyze(arg)


    one = time.process_time() 
    print(" | Time elapsed: ", one , "sec")

    # for File in FileNames:
    #     ImportDataFromHDF5(File, channels)
    # for amp in ch1.Amp: 
    #     test = np.max(amp) 
    #     # print(np.max(amp))
    # for amp in ch2.Amp: 
    #     test = np.max(amp) 
    #     # print(np.max(amp))



    two = time.process_time() 
    print(" | Time elapsed: ", two-one , "sec")

    for ch in channels: 
        for Files in FileNames: 
            test = Yield(Files)
            np.mean(Yield(Files)[0])
            for i in Yield(Files): 
                test = np.max(i[0])

            
    three = time.process_time() 
    print(" | Time elapsed: ", three-two , "sec")
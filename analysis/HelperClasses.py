import datetime
import os
import glob
import time
import sys

def ProgressBar(it, counts):
    width = counts
    if(it==1):
        sys.stdout.write(" | [%s]" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width+1))
        sys.stdout.write("%")
    else:
        sys.stdout.write("%")
        sys.stdout.flush()
    if(it==counts):
        sys.stdout.write("\n")

def ElapsedTime():
    print " | Time elapsed:          ", time.clock() , "sec"

class MetaData:
    def __init__(self, InputPath):
        self.InputPath = InputPath
        # self.DateOfDataTaking = self.InputPath.split('/')[-3]
        # self.ExactTimeOfDataTaking = datetime.datetime.strptime(time.ctime(os.path.getctime(InputPath)), "%a %b %d %H:%M:%S %Y").strftime('%Y%m%d%H%M%S')
        self.DateOfProcessing = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.ProjectPath = '/project/fas/david_moore/aj487/LXeData/'
        # self.FileNames = glob.glob(self.InputPath+"*.txt")
        # self.NrOfFiles = len(self.FileNames)
        # self.OutputPath = self.CreateDirectory()
        # self.DataName = self.FileNames[0].split('/')[-2]


    def CreateDirectory(self):
        Directory = self.ProjectPath + self.DateOfDataTaking + '/'
        if not os.path.exists(Directory):
            os.makedirs(Directory)
        return Directory

import datetime
import os
import glob
import time

class MetaData:
    def __init__(self, InputPath):
        self.InputPath = InputPath
        self.DateOfDataTaking = self.InputPath.split('/')[-3]
        self.ExactTimeOfDataTaking = datetime.datetime.strptime(time.ctime(os.path.getctime(InputPath)), "%a %b %d %H:%M:%S %Y").strftime('%Y%m%d%H%M%S')
        self.DateOfProcessing = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        self.ProjectPath = '/project/fas/david_moore/aj487/LXeData/'
        self.FileNames = glob.glob(self.InputPath+"*.txt")
        self.NrOfFiles = len(self.FileNames)
        self.OutputPath = self.CreateDirectory()
        self.DataName = self.FileNames[0].split('/')[-2]


    def CreateDirectory(self):
        Directory = self.ProjectPath + self.DateOfDataTaking + '/'
        if not os.path.exists(Directory):
            os.makedirs(Directory)
        return Directory

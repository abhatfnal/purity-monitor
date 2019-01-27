import datetime
import os

class MetaData:
    def __init__(self, InputPath):
        self.InputPath = InputPath
        self.DateOfDataTaking = self.InputPath.split('/')[-3]
        self.DateOfProcessing = datetime.date.today().strftime("%Y%m%d")
        self.ProjectPath = '/project/fas/david_moore/aj487/LXeData/'
        self.FileNames = glob.glob(options.filepath+"*.txt")
        self.NrOfFiles = len(self.FileNames)
        self.OutputPath = self.CreateDirectory()

    def CreateDirectory(self):
        dirname = self.ProjectPath + self.DateOfDataTaking + '/'
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        return dirname 

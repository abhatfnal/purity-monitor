import glob
from WaveformClass import WFM


def GetFileCreationTime(file):
    timestamp = os.path.getmtime(file)
    modificationTime = datetime.datetime.fromtimestamp(timestamp)
    return modificationTime

def ReadFilesInDirectories(filepath):
    directories = sorted(os.listdir(filepath))
    files = []
    NrFiles = 0
    for dir in directories:
        NewDir = filepath+"/"+dir+"/"
        print " | Path:\t", NewDir
        NrFiles2, Files = ReadFilesInDirectory(NewDir)
        files.append(Files)
        NrFiles += NrFiles2
    files = [val for sublist in files for val in sublist]
    return NrFiles, files

def ReadFilesInDirectory(filepath):
    files = glob.glob(filepath+"*.txt")
    files2 = sorted(files)
    NrFiles = len(files)
    print " | Number of files:\t", NrFiles
    return NrFiles, files2

def ReadData(channels, files):
    print " | Reading in data files..."
    for i,line in enumerate(files):
        # ProgressBar(i+1, len(files))
        with open(line) as wfm:
            for j,data in enumerate(wfm):
                columns = data.split()
                for k, ch in enumerate(channels):
                    if(j==0):
                        ch.CreateTime.append(GetFileCreationTime(line))
                    if(i==0):
                        ch.Time.append(float(columns[0])*ch.TScale)
                    ch.Amp[i].append(float(columns[ch.ID])*ch.VScale)

def TxtInput(options, channels):
    if(options.txt):
        #Read in files in a given directory. Filename is the file that contains the names of all the data files. num_lines is the number of such data files.
        if(options.dir):
            num_files, files = ReadFilesInDirectories(options.filepath)
        else:
            num_files, files = ReadFilesInDirectory(options.filepath)
        #Get all waveforms in data files and save them in lists. The content of each waveform is saved in chX.Amp[i], with X the channel number and i the waveform number. The time is saved only once in chX.Time
        ReadData([ch1, ch2], files)
        if(options.filepath2 != ""):
            num_files2, files2 = ReadFilesInDirectory(options.filepath2)
            ch3 = WFM(num_files2, options.filepath2, 2)
            ReadData([ch3], files2)
    else:
        #New Method of getting data from HDF5 files.
        ImportDataFromHDF5(options.filepath, [ch2])

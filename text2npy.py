import time
import h5py
from optparse import OptionParser
from wv_class import WFM
from helper_classes import MetaData
from features import ProgressBar

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
(options, args) = parser.parse_args()

def ReadData(channels, meta):
    print " | Number of files...        ", meta.NrOfFiles
    print " | Reading in data files..."
    for i,line in enumerate(meta.FileNames):
        ProgressBar(i+1, meta.NrOfFiles)
        wfm = open(line)
        for j,data in enumerate(wfm):
            columns = data.split()
            for k, ch in enumerate(channels):
                if(i==0):
                    ch.Time.append(float(columns[0])*ch.TScale)
                ch.Amp[i].append(float(columns[ch.ID])*ch.VScale)
        wfm.close()

def WriteToHDF5(channels, meta):
    print " | Writing data to HDF5 file..."
    with h5py.File(meta.OutputPath + meta.DataName + '.h5', 'w') as hf:
        g = hf.create_group('Default')
        g.attrs['NrOfFiles'] = meta.NrOfFiles
        g.attrs['DateOfProcessing'] = meta.DateOfProcessing
        g.attrs['DateOfDataTaking'] = meta.DateOfDataTaking
        g.attrs['ExactTimeOfDataTaking'] = meta.ExactTimeOfDataTaking
        g.attrs['DataName'] = meta.DataName
        g.attrs['InputPath'] = meta.InputPath
        g.create_dataset('Time', data=channels[0].Time)
        for i in range(meta.NrOfFiles):
            g.create_dataset('Waveforms/Waveform%d' % i, data=channels[1].Amp[i])
        for ch in channels:
            ch.GetSampling(state=True)
            ch.GetAverageWfm(state=True)
        g.create_dataset('Waveforms/Average', data=channels[1].MeanAmp)
        g.create_dataset('Trigger', data=channels[0].MeanAmp)


if __name__ == '__main__':
    Meta = MetaData(options.filepath)
    ch1 = WFM(Meta.NrOfFiles, 1, "Trigger", 1)
    ch2 = WFM(Meta.NrOfFiles, 1, "Signal", 2)

    ReadData([ch1,ch2], Meta)
    WriteToHDF5([ch1,ch2], Meta)

    print " | Time elapsed: ", time.clock() , "sec"

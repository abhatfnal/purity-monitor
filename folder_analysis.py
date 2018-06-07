from optparse import OptionParser
import time, datetime, sys, os, glob, struct

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
(options, args) = parser.parse_args()

def ReadFilesInDirectory():
    filename = "dir.dat"
    os.system("ls "+options.filepath+" > "+filename)
    num_files = sum(1 for line in open(filename))
    print " | Number of directories:          ", num_files
    return num_files, filename

def RunAnalysis(filename):
    file = open(filename, "r")
    for i,line in enumerate(file):
        print " | Analizing...", os.path.join(line[:-1], '')
        # print ("python read.py -f "+options.filepath+os.path.join(line, ''))
        os.system("python read.py -f "+options.filepath+os.path.join(line[:-1], ''))

if __name__ == '__main__':
    num_files, filename = ReadFilesInDirectory()
    RunAnalysis(filename)

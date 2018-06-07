from optparse import OptionParser
import os

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--file", action="store", type="string", dest="filepath", help="Specify path/file to be read in.")
parser.add_option("-m", "--mail", action="store", type="string", dest="mail", help="Run each file as batch job.")
(options, args) = parser.parse_args()

def ReadFilesInDirectory():
    filename = "dir.dat"
    os.system("ls "+options.filepath+" > "+filename)
    filenumber = sum(1 for line in open(filename))
    print " | Number of directories:          ", filenumber
    return filenumber, filename

def RunAnalysis(filename, filenumber):
    file = open(filename, "r")
    for i,line in enumerate(file):
        print " | Analizing...", os.path.join(line[:-1], '')
        os.system("python read.py -f "+options.filepath+os.path.join(line[:-1], ''))

def RunCluster(filename, filenumber):
    file = open(filename, "r")
    for i,line in enumerate(file):
        f = open("charge.sh","w")
        f.write("#!/bin/bash\n")
        if(i==0):
            f.write("#SBATCH --mail-type=BEGIN\n")
            f.write("#SBATCH --mail-user="+options.mail+"\n")
        if(i==filenumber-1):
            f.write("#SBATCH --mail-type=END\n")
            f.write("#SBATCH --mail-user="+options.mail+"\n")
        f.write("#SBATCH --partition=day\n")
        f.write("#SBATCH --job-name=cluster\n")
        f.write("#SBATCH --ntasks=1 --nodes=1\n")
        f.write("#SBATCH --mem-per-cpu=4000\n")
        f.write("#SBATCH --time=0:01:00\n\n")
        f.write("python read.py -f "+options.filepath+os.path.join(line[:-1], ''))
        f.close()
        os.system("sbatch charge.sh")
        if(i==filenumber-1):
            break

if __name__ == '__main__':
    filenumber, filename = ReadFilesInDirectory()
    if not options.mail:
        RunAnalysis(filename, filenumber)
    else:
        RunCluster(filename, filenumber)

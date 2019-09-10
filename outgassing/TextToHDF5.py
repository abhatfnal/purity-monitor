import h5py, argparse, datetime, os, time, sys, decimal

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, action="store",  dest="input", nargs="*")
parser.add_argument("-o", "--output", type=str, action="store",  dest="output")
arg = parser.parse_args()

def ProgressBar(count, total, phrase):
    rows, columns = os.popen('stty size', 'r').read().split()
    barLength = 30 # Modify this to change the length of the progress bar
    status = '%d/%d' % (count+1, total)
    progress = count/float(total)
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if count == total-1:
        progress = 1
        status = status+"           \n"
    block = int(round(barLength*progress))
    text = '\r%-*s [%s] %-*s %s' % (25,phrase,"#"*block+"-"*(barLength-block),8,'%.2f%%' % decimal.Decimal(progress*100),status)
    sys.stdout.write(text)
    sys.stdout.flush()

def ChooseFilesToAnalyze(arg):
    files = []
    if (arg.input != None):
        if ('*' in arg.input):
            files.append(glob.glob(dir+arg.input))
        else:
            files.append(arg.input) 	
    files = [val for sublist in files for val in sublist]
    return files

if __name__ == '__main__':
    Files = ChooseFilesToAnalyze(arg)
    Outputpath = '/project/fas/david_moore/zl423/PurityData/Outgassing/'
    h1 = h5py.File(Outputpath+arg.output, 'w')
    g1 = h1.create_group('Pressure')
    for jj, File in enumerate(Files): 
        ProgressBar(jj, len(Files), 'Converting txt files')
        f1 = open(File, "r")
        mass = [] 
        pressure = [] 
        for ii, line in enumerate(f1): 
            if ii < 6: 
                pass 
            if ii > 5 and ii < 20: 
                data = line.split(',')
                h1.attrs[u'%s' % data[0]] = data[1]
            if ii > 21:
                data = line.split(',')
                mass.append(float(data[0]))
                pressure.append(float(data[1]))
        path, filename = os.path.split(File)
        if('AM' in filename or 'PM' in filename): 
            dt = datetime.datetime.strptime(filename, '%b_%d_%Y__%I-%M-%S_%p.txt')
            do = datetime.datetime.strftime(dt, "%Y%m%d%H%M%S")
        g1.create_dataset(do,data=pressure)
    h1.create_dataset('Mass',data=mass)
    h1.close()
    print("Time elapsed: ", time.process_time() , "sec")
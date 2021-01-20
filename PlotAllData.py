import os
import time
import glob 

files = []
path = '/project/fas/david_moore/zl423/PurityData/StandStatus/'
files = glob.glob(path+"*.h5")
for name in files: 
    os.system("python StandMonitor.py -f %s" % (name))
    print time.clock()

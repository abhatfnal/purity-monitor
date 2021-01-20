import os
import sys
import time

path = '/project/fas/david_moore/zl423/PurityData/StandStatus/'
name = '20191206.h5'

file = path+name

keep_time = os.stat(file).st_mtime

while True:
    time.sleep(30)
    second_time = os.stat(file).st_mtime
    if(second_time == keep_time):
        continue
    else:
        try:
            os.system("python StandMonitor.py -f %s" % (file))
        except:
            continue
    keep_time = second_time

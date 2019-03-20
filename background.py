import os
import sys
import time

path = '/project/fas/david_moore/zl423/PurityData/StandStatus/'
name = '3_19_2019'
date = '20190319'

file = path+name

keep_time = os.stat(file).st_mtime

while True:
    time.sleep(30)
    second_time = os.stat(file).st_mtime
    if(second_time == keep_time):
        continue
    else:
        os.system("python monitor.py -f %s -d %s" % (file, date))
    keep_time = second_time

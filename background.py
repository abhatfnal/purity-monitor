import time
import os


file = "$PDATA/StandStatus/1_14_2019"
date = "20190114"
while True:
    os.system("python monitor.py -f %s -d %s" % (file, date))
    time.sleep(60)
    print 'next round'

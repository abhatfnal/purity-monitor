# -*- coding: utf-8 -*-
#!/usr/bin/python

import sys
import redpitaya_scpi as scpi
import matplotlib.pyplot as plt
import struct

num_traces = 100 #number of traces to record
dec_fac = 1 # decimation factor (see https://redpitaya.readthedocs.io/en/latest/appsFeatures/examples/acqRF-samp-and-dec.html)
trig_lev_mv = 25 # trigger level in mV

rp_s = scpi.scpi('169.254.39.61')

rp_s.tx_txt('ACQ:DATA:FORMAT BIN')
rp_s.tx_txt('ACQ:DATA:UNITS RAW')
rp_s.tx_txt('ACQ:DEC %d'%dec_fac)
rp_s.tx_txt('ACQ:AVG OFF')

rp_s.tx_txt('ACQ:TRIG:LEV %d mV'%trig_lev_mv)
rp_s.tx_txt('ACQ:TRIG CH1_PE')
rp_s.tx_txt('ACQ:TRIG:DLY 8000')
rp_s.tx_txt('ACQ:START')

i=0
data_vec = []
while 1:
    rp_s.tx_txt('ACQ:TRIG:STAT?')
    if rp_s.rx_txt() == 'TD':
        i += 1
        print(i)
        rp_s.tx_txt('ACQ:SOUR1:DATA?')
        buff_byte = rp_s.rx_arb()
        data_vec.append(buff_byte)
        rp_s.tx_txt('ACQ:TRIG CH1_PE;START')
        if(i>100): break


plt.figure()
for pidx,d in enumerate(data_vec):
    buff = [struct.unpack('!h',bytearray(d[i:i+2]))[0] for i in range(0, len(d), 2)]
    if(pidx < 10):
        plt.plot(buff, label=str(pidx))
plt.ylabel('Voltage')

plt.show()
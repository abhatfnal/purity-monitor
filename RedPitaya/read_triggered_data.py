# -*- coding: utf-8 -*-
#!/usr/bin/python

import redpitaya_scpi as scpi
import matplotlib.pyplot as plt
import struct
import numpy as np
import h5py

rp_s = scpi.scpi('169.254.84.97')

def acquire_file(num_traces=1000, dec_fac=2, trig_lev_mv=25):
    
    #num_traces = 100 #number of traces to record
    #dec_fac = 1 # decimation factor (see https://redpitaya.readthedocs.io/en/latest/appsFeatures/examples/acqRF-samp-and-dec.html)
    #trig_lev_mv = 25 # trigger level in mV
    
    rp_s.tx_txt('ACQ:DATA:FORMAT BIN')
    rp_s.tx_txt('ACQ:DATA:UNITS RAW')
    rp_s.tx_txt('ACQ:DEC %d'%dec_fac)
    rp_s.tx_txt('ACQ:AVG OFF')
    rp_s.tx_txt('ACQ:TRIG:LEV %d mV'%trig_lev_mv)
    rp_s.tx_txt('ACQ:TRIG CH1_PE')
    rp_s.tx_txt('ACQ:TRIG:DLY 8000')
    rp_s.tx_txt('ACQ:START')
    
    i=0
    data_vec_ch1 = []
    data_vec_ch2 = []
    while 1:
        rp_s.tx_txt('ACQ:TRIG:STAT?')
        if rp_s.rx_txt() == 'TD':
            i += 1
            if(i%100==0): print(i)
            rp_s.tx_txt('ACQ:SOUR1:DATA:OLD:N? 4096')
            buff_byte = rp_s.rx_arb()
            data_vec_ch1.append(buff_byte)
            rp_s.tx_txt('ACQ:SOUR2:DATA:OLD:N? 4096')
            buff_byte = rp_s.rx_arb()
            data_vec_ch2.append(buff_byte)
            rp_s.tx_txt('ACQ:TRIG CH1_PE;START')
            if(i>=num_traces): break
    
    ch1_dat = []
    for pidx,d in enumerate(data_vec_ch1):
        buff = [struct.unpack('!h',bytearray(d[i:i+2]))[0] for i in range(0, len(d), 2)]
        ch1_dat.append(buff)
    
    ch2_dat = []
    for pidx,d in enumerate(data_vec_ch2):
        buff = [struct.unpack('!h',bytearray(d[i:i+2]))[0] for i in range(0, len(d), 2)]
        ch2_dat.append(buff)
    
    return np.array(ch1_dat), np.array(ch2_dat)

cd1, cd2 = acquire_file(num_traces=100)

channel1=[]
channel2=[]
for x,y in zip(cd1.T,cd2.T):
    channel1.append(x*2/16384)
    channel2.append(y*2/16384)

f = h5py.File("mytestfile.h5", "w")
f.create_dataset('Ch1',data=channel2)
f.create_dataset('Ch2',data=channel2)
f.close()
plt.figure()
plt.plot(channel1,color='k')
plt.plot(channel2,color='r')
# plt.legend(loc='upper right')
plt.show()

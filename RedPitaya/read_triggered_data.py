# -*- coding: utf-8 -*-
#!/usr/bin/python

import redpitaya_scpi as scpi
import matplotlib.pyplot as plt
import struct
import numpy as np
import h5py
from datetime import datetime as dt

start_time=dt.now()

# Sampling rate of redpitaya's ADC is always 125MS/s per channel, so it takes 1/(125*10^6)=8 ns to take 1 sample
#Max Buffer length=16384
#To take a decimated sample it takes dec_fac*8ns
#To fill the buffer with 16384 samples, the time is t=16384*dec_fac*8ns

def acquire_file(rp_s, decimation, buffer_size, trg_delay, trig_channel, num_traces=10, dec_fac=1, trig_lev_mv=100):
    
    #num_traces = 100 #number of traces to record
    #dec_fac = 1 # decimation factor (see https://redpitaya.readthedocs.io/en/latest/appsFeatures/examples/acqRF-samp-and-dec.html)
    #trig_lev_mv = 25 # trigger level in mV
    
    rp_s.tx_txt('ACQ:DATA:FORMAT BIN')
    rp_s.tx_txt('ACQ:DATA:UNITS RAW')
    rp_s.tx_txt('ACQ:DEC %d'%dec_fac)
    rp_s.tx_txt('ACQ:AVG OFF')

    rp_s.tx_txt('ACQ:START')

    rp_s.tx_txt('ACQ:TRIG CH%d_PE'%trig_channel)
    # rp_s.tx_txt('ACQ:TRIG NOW')
    rp_s.tx_txt('ACQ:TRIG:LEV %d mV'%trig_lev_mv)
    rp_s.tx_txt('ACQ:TRIG:DLY %d'%trg_delay)
    rp_s.tx_txt('ACQ:TRIG:DLY?')
    # rp_s.tx_txt('ACQ:TPOS?')
    # print('-'*10)
    # print('TPOS: ',rp_s.rx_txt())

    
    i=0
    data_vec_ch1 = []
    data_vec_ch2 = []
    while 1:
        rp_s.tx_txt('ACQ:TRIG:STAT?')
        check = rp_s.rx_txt()
        # print(i, check )
        if check == 'TD':
            i += 1
            
            if(i%100==0): print(i)
            # rp_s.tx_txt('ACQ:SOUR1:DATA:OLD:N? %d'%buffer_size)
            # buff_byte = rp_s.rx_arb()
            # data_vec_ch1.append(buff_byte)
            rp_s.tx_txt('ACQ:SOUR2:DATA:OLD:N? %d'%buffer_size)
            # rp_s.tx_txt('ACQ:SOUR2:DATA?')
            buff_byte = rp_s.rx_arb()
            data_vec_ch2.append(buff_byte)
            # rp_s.tx_txt('ACQ:TRIG CH%d_PE;START'%trig_channel)
            # rp_s.tx_txt('ACQ:START')
           
            if(i>=num_traces): break
    
    # ch1_dat = []
    # for pidx,d in enumerate(data_vec_ch1):
    #     buff = [struct.unpack('!h',bytearray(d[i:i+2]))[0] for i in range(0, len(d), 2)]
    #     ch1_dat.append(buff)
    
    ch2_dat = []
    for pidx,d in enumerate(data_vec_ch2):
        buff = [struct.unpack('!h',bytearray(d[i:i+2]))[0] for i in range(0, len(d), 2)]
        ch2_dat.append(buff)
        # print(len(buff))

    return np.array(ch2_dat), np.array(ch2_dat)


rp_s = scpi.scpi('172.28.175.57')
decimation=2**2 #should be an integer which is 2^n like 1,2,4,8,16,32 etc.
max_buffer_size=16384
buffer_size=1024*2
# buffer_size=max_buffer_size
trg_delay=8196-(buffer_size/2.0)
# trg_delay=0
time_waveform =np.linspace(0,buffer_size*8E-9*decimation,buffer_size)
trig_channel=2 #can be ony 1 or 2


# plt.figure()
# channel1=[]
channel2=[]
num_wvf=100
for i in range(num_wvf):
    # print(i)
    cd1,cd2=acquire_file(rp_s, decimation, buffer_size, trg_delay, trig_channel, num_traces=1, dec_fac=decimation, trig_lev_mv=50)
    # print(cd2[0])
    channel2.append(np.array(cd2)[0]*1000*2/max_buffer_size)
    # plt.plot(time_waveform*1E6, np.array(cd2)[0]*1000*2/max_buffer_size ,lw=0.5, label='{}'.format(i))

# plt.xlabel('Time [$\mu$s]')
# plt.ylabel('Amplitude [mV]')
# # plt.legend(loc='upper left')
# plt.show()






# for x,y in zip(cd1.T,cd2.T):
#     channel1.append(1000*x*2/max_buffer_size) #For LV, Maximum Voltage is +-1V
#     channel2.append(1000*y*2/max_buffer_size)
#     # channel1.append(x*40/max_buffer_size) #For HV, Maximum Voltage is +-20V
#     # channel2.append(y*40/max_buffer_size)
#     # print(1000*y*2/max_buffer_size)
f = h5py.File("mytestfile.h5", "w")
# f.create_dataset('Ch1',data=channel1)
f.create_dataset('Ch2',data=channel2)
f.create_dataset('Time',data=time_waveform)
f.close()


# plt.figure()
# # plt.plot(time*1E6,channel1,color='k',lw=0.5)
# # plt.plot(time_waveform*1E6,np.array(channel2)[:,0],lw=0.5)
# # plt.plot(time_waveform*1E6+10,np.array(channel2)[:,1],lw=0.5)
# # plt.plot(time_waveform*1E6+20,np.array(channel2)[:,2],lw=0.5)
# for i in range(num_wvf):
#     plt.plot(time_waveform*1E6,np.array(channel2)[:,i],lw=0.5, label='{}'.format(i))

# plt.xlabel('Time [$\mu$s]')
# plt.ylabel('Amplitude [mV]')
# plt.legend(loc='upper left')
# plt.show()
end_time=dt.now()


print('Duration: {}'.format(end_time-start_time))

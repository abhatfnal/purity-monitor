

N = num_points
T = 1.0 / 5000000.0
x = np.linspace(0.0, N*T, N)
yf = fft(ch1)
yf2 = fft(ch2)
w = blackman(N)
ywf = fft(ch1*w)
ywf2 = fft(ch2*w)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)


fs = 5000000.0
lowcut = 0
highcut = 100E3
T = 0.002
nsamples = T * fs
t = np.linspace(-T, T, nsamples, endpoint=False)

_, y = GetFFTWaveform(data=ch1,lowcut=0,highcut=200E3,sampling=5000000.0,points=num_points,order=6)
_, y2 = GetFFTWaveform(data=ch2,lowcut=0,highcut=200E3,sampling=5000000.0,points=num_points,order=6)
yf3 = fft(y)
yf4= fft(y2)


start = -100
end = 15

start1 = 100
end1 = 900

start2 = time[y.index(min(y))]
end2 = 900

start3 = time[y2.index(max(y2))]
end3 = 750



# popt, pcov = curve_fit(exponential_func, np.asarray(time[time.index(start):time.index(end)]), y[time.index(start):time.index(end)], p0=(-1, 1e-6, 0, 1))

# popt1, pcov1 = curve_fit(exponential_func, np.asarray(time[time.index(start1):time.index(end1)]), y[time.index(start1):time.index(end1)], p0=(-.1, 1e-4, 1, 1))
#
# popt2, pcov2 = curve_fit(double_exponential_func, np.asarray(time[time.index(start2):time.index(end2)]), y[time.index(start2):time.index(end2)], p0=(-1.0, 0.1, -1.0, 0.1, 1.0))
#
# popt3, pcov3 = curve_fit(exponential_func, np.asarray(time[time.index(start3):time.index(end3)]), y2[time.index(start3):time.index(end3)], p0=(1, 1e-4, 0, 1))


# popt, pcov = curve_fit(triple_exponenial_func, np.asarray(time[time.index(start):time.index(end)]), y[time.index(start):time.index(end)], p0=(1, 1e-6, 1, 1, 1e-6,1, 1e-6, 1,1,1))
#
#
#



# yy = exponential_func(np.asarray(time), *popt)
# yy1 = exponential_func(np.asarray(time), *popt1)
#
# yy2 = double_exponential_func(np.asarray(time), *popt2)
# yy3 = exponential_func(np.asarray(time), *popt3)


# print "Parameters:  ", popt
# print "Parameters2:  ", popt2
# print "Cathode Value Double at 0:  ", double_exponential_func(time[y.index(min(y))+1],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4])
# idx = (np.abs(y-value)).argmin()
# print "Cathode Value at 0:\t", exponential_func(time[y.index(min(y))],popt1[0],popt1[1],popt1[2],popt1[3]), '\t', min(y), '\t', time[y.index(min(y))]
# print "Anode Value at 0:\t", exponential_func(time[y2.index(max(y2))],popt3[0],popt3[1],popt3[2],popt3[3]), max(y2), '\t', time[y2.index(max(y2))]



# sum = [x+z for x,z in zip(y,y2)]

fig = plt.figure(figsize=(12,7))
ax = fig.gca()
ax.grid()
plt.xlabel("Time [us]")
plt.ylabel("Voltage [mV]")
plt.xlim(-200,200)
plt.ylim(-100,100)
plt.plot(time, y)
plt.plot(time, y2)
# plt.plot(time, yy)
# plt.plot(time, yy3)
# plt.plot(time, yy1)
# plt.plot(time, yy2)
plt.legend(['Cathode', 'Anode','Falling Anode','Falling Cathode', 'Falling Cathode Double'])
plt.show()




fig = plt.figure(figsize=(12,7))

plt.subplot(421)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, ch1, '-b')

plt.subplot(423)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, ch2, '-r')

plt.subplot(422)
plt.grid(True)
plt.title("Frequency band pass filter between 0 and 100 kHz")
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2], 2.0/N * np.abs(yf[1:N//2]), '-b')

plt.subplot(424)
plt.grid(True)
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf2[1:N//2]), '-r')

plt.subplot(425)
plt.grid(True)
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3, y, '-b')

plt.subplot(427)
plt.grid(True)
plt.xlabel("Time [ms]")
plt.ylabel("Amplitude [mV]")
plt.plot(t*1E3,y2, '-r')

plt.subplot(426)
plt.grid(True)
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf3[1:N//2]), '-b')

plt.subplot(428)
plt.grid(True)
plt.xlabel("Frequency [kHz]")
plt.ylabel("Amplitude")
plt.semilogy(xf[1:N//2]/1E3, 2.0/N * np.abs(yf4[1:N//2]), '-r')

plt.show()

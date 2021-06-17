import numpy as np
import matplotlib.pyplot as plt
import scipy

def GetExp(x,a,b):
    return a*np.exp(-b*x)

def OptimalFilter(Time, Data, Pol):
    Time = np.array(Time)
    SamplePoints = Time.size
    SampleSpace = 1/(1E6/(Time[-1]-Time[-2]))
    TimeFFT = np.fft.fftfreq(SamplePoints, SampleSpace)[0:SamplePoints//2]
    print SamplePoints, SampleSpace, SamplePoints*SampleSpace
    
    Data = [x*Pol for x in Data]
    Template1 = np.mean(Data, axis=0)
    Template = Template1/np.max(Template1)
    TemplateMaxBin = np.where(Template==np.max(Template))[0][0]
    binof = 2000
    binoff = 10000
    First = TemplateMaxBin+binof
    Second = TemplateMaxBin+binof+binoff
    x = range(First, Second)
    x = np.linspace(First, Second, Second-First)
    print x 
    
    amp = []
    fall = [] 
    for data in Data:
        init_vals = [np.max(data),0.001]
        p,_ = scipy.optimize.curve_fit(GetExp, Time[First:Second], data[First:Second], p0=init_vals)
        # p,_ = scipy.optimize.curve_fit(lambda x,a,b: GetExp(x,np.max(data),b), x, data[TemplateMaxBin+100:] )
        print p[0],p[1]
        plt.plot(Time, data, color='black')
        # plt.plot(Time[First:Second], GetExp(Time[First:Second], *p), color='red')
        # plt.ion()
        # plt.draw()
        amp.append(p[0])
        fall.append(p[1])
    
    # plt.scatter(amp,fall)
    plt.plot(Time, Template1, color='red', label='Template')
    plt.show()
    quit()
    
    plt.xlabel('Time [us]')
    plt.ylabel('Amplitude [mV]')
    plt.legend(loc='upper left')
    plt.xlim(-1000,1000)
    plt.grid()
    plt.show()
    quit() 
    TemplateMaxBin = np.where(Template==np.max(Template))[0][0]
    TDiff = [np.where(x==np.max(x))[0][0] - TemplateMaxBin for x in Data]
    # TDiff = [np.where(x==np.max(x))[0][0] for x in Data]
    # TDiff = [Time[np.where(x==np.max(x))[0][0]] - Time[TemplateMaxBin] for x in Data]

    NoiseWaveforms = GetNoiseWaveforms(Data)
    NoiseFFT = GetFFT(NoiseWaveforms, SamplePoints)
    PowerSpectrumDensity = GetPowerSpectrumDensity(NoiseWaveforms, 1/SampleSpace,SamplePoints)

    TemplateFFT = GetFFT([Template], SamplePoints)[0]
    TemplateFFTConj = np.conj(TemplateFFT)

    DataFFT = GetFFT(Data, SamplePoints)
    AmplitudeEstimator = GetAmplitudeEstimator(Data, DataFFT, TemplateFFT, TemplateFFTConj, PowerSpectrumDensity, TDiff)
    
    # OptFilt = np.fft.irfft(TemplateFFTConj*DataFFT[0]/PowerSpectrumDensity)
    # OptFilt = np.fft.irfft(TemplateFFTConj/PowerSpectrumDensity)
    # OptFilt = TemplateFFTConj*DataFFT[0]/PowerSpectrumDensity

    amplitude = []
    for ii, (x,y) in enumerate(zip(Data,NoiseWaveforms)):
        OptFilt = TemplateFFTConj*DataFFT[ii]/PowerSpectrumDensity
        fig = plt.figure(figsize=(16,9))
        plt.subplot(121)
        print np.shape(TimeFFT), np.shape(DataFFT[ii])
        plt.plot(DataFFT[ii], label='Data FFT', color='blue')
        plt.plot(TemplateFFTConj, label='Template FFT Conj', color='green')
        plt.plot(TemplateFFT, label='Template FFT ', color='pink')
        plt.plot(OptFilt, label='Optimal Filter Data FFT', color='red')
        plt.plot(NoiseFFT[ii], label='Noise FFT', color='black')
        plt.plot(TemplateFFTConj*DataFFT[ii], label='Template FFT Conj * Data FFT', color='grey')
        plt.plot(PowerSpectrumDensity, label='Noise PSD', color='orange')
        plt.xlabel('index')
        plt.ylabel('Intensity')
        plt.xlim(1,len(PowerSpectrumDensity))
        # plt.ylim(1E-10,1E7)
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.subplot(122)
        plt.xlabel('Time [us]')
        plt.ylabel('Amplitude [mV]')
        plt.plot(Time, x, label='Data', color='blue')
        plt.plot(Time, Template1, label='Template', color='green')
        plt.plot(Time, Template*AmplitudeEstimator[ii], label='New', color='green')
        new = np.fft.irfft(OptFilt)
        print new 
        new = new/(np.max(new)/np.max(Data[ii]))
        print new 
        plt.plot(Time, new, label='Optimally Filtered Data', color='red')
        plt.plot(Time, y, label='Noise', color='black')
        # plt.plot(Time, x-y, label='Subtracted', color='yellow')
        plt.legend(loc='upper left')
        plt.xlim(-1000,1000)
        plt.show()

        amplitude.append(np.max(x - y))
    
    print np.array(amplitude)
    quit()
    maxorigin = np.array([np.max(x) for x in Data])
    print np.array(amplitude) - maxorigin
    bins = np.linspace(15,25,100)
    plt.hist(amplitude, bins, alpha = 0.5)
    plt.hist(maxorigin, bins, alpha = 0.5)
    plt.show()
    quit()
    return AmplitudeEstimator

def GetPowerSpectrumDensity(Noise, Sampling, SamplePoints):
    NoiseFFT = [np.fft.rfft(x) for x in Noise]

    # PowerSpectrumDensity = []
    # for ii,x in enumerate(Noise): 
    #     f, Pxx = scipy.signal.periodogram(x, Sampling)
    #     PowerSpectrumDensity.append(Pxx)
    # PowerSpectrumDensity = np.array(PowerSpectrumDensity)
    # PowerSpectrumDensity = np.mean(PowerSpectrumDensity, axis=0)

    PowerSpectrumDensity = [np.abs(x)**2 for x in NoiseFFT]
    PowerSpectrumDensity = np.array(PowerSpectrumDensity)
    PowerSpectrumDensity = np.mean(PowerSpectrumDensity, axis=0)
    return PowerSpectrumDensity

def GetNoiseWaveforms(Data):
    Template = np.mean(Data, axis=0)
    Noise = np.array([x-Template/(np.max(Template)/np.max(x)) for x in Data])
    return Noise

def GetFunc(x, A, Template): 
    func = A * Template
    return func 

def GetAmplitudeEstimator(Data, DataFFT, TemplateFFT, TemplateFFTConj, PowerSpectrumDensity, TDiff):
    Amplitudes = []
    NrPoints = np.shape(Data)[1]/2
    print NrPoints
    j = complex(0, 1)
    best = []
    omega = np.linspace(-2*np.pi/NrPoints, 2*np.pi/NrPoints, NrPoints)
    for ii, datafft in enumerate(DataFFT):   
        new = [] 
        for t0 in range(0,1):
            Numerator = np.sum( np.exp(j*omega*t0) * TemplateFFTConj[1:] * datafft[1:] / PowerSpectrumDensity[1:] ).real
            # Numerator = np.sum(TemplateFFTConj[1:] * datafft / PowerSpectrumDensity ).real
            # print np.sum(TemplateFFTConj * datafft / PowerSpectrumDensity )
            # print np.abs(np.sum(TemplateFFTConj * datafft / PowerSpectrumDensity ).real)/np.abs(np.sum(TemplateFFTConj * datafft / PowerSpectrumDensity ).imag)
            Denominator = np.sum((np.abs(TemplateFFT[1:])**2)/PowerSpectrumDensity[1:])
            AmplitudeEstimator = Numerator/Denominator
            new.append(AmplitudeEstimator)
        # plt.figure()
        # plt.plot(range(-NrPoints,NrPoints,NrPoints/100), new )
   
        p,_ = scipy.optimize.curve_fit(lambda x,a: GetFunc(x,a, np.mean(Data, axis=0)[60000:100000]/np.max(np.mean(Data, axis=0))), range(-NrPoints,NrPoints,NrPoints/100), Data[ii][60000:100000], )
        # plt.figure()
        # plt.plot(Data[ii], color='black')
        # plt.plot(np.mean(Data, axis=0)/np.max(np.mean(Data, axis=0))*p,color='green')
        # plt.plot(np.mean(Data, axis=0)/np.max(np.mean(Data, axis=0))*np.max(new),color='red')
        # plt.show()
        best.append(p)
        print AmplitudeEstimator
        Amplitudes.append(AmplitudeEstimator)
    

    Amplitudes = np.array(Amplitudes)
    maxorigin = np.array([np.max(x) for x in Data])
    plt.hist(best, 100)
    plt.hist(maxorigin, 100, color='red')
    # plt.hist(maxorigin*Amplitudes, 100, color='green')
    # plt.show()
    return Amplitudes

def GetFFT(Data, SamplePoints): 
    # DataFFT = [ 2.0/SamplePoints * np.abs(np.fft.rfft(x))[0:SamplePoints//2] for x in Data]
    DataFFT = [np.fft.rfft(x) for x in Data]
    return np.array(DataFFT)

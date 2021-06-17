import numpy as np
import matplotlib.pyplot as plt

def MatchedFilter(Time, Data): 
    fft_cut = 4000

    Template = np.mean(Data, axis=0)

    new = vary_pulse_time(Data[0], Time, Template)
    print(new)
    # plt.plot(Data[0])
    plt.plot(new)
    plt.yscale('log')
    plt.show()
    quit() 


    fft_template_fit = np.fft.rfft(Template)    
    fft_template_fit[fft_cut:].fill(0.0)
    fft_template_fit_conj = np.conj(fft_template_fit)
    ifft_template_fit = np.fft.irfft(fft_template_fit).real
    # plot_fft(fft=[np.abs(fft_template_fit)**2])
    # plot_wfs(wfs=[Template, ifft_template_fit], label=['fit', 'ifft_fit'])
    # plt.show()
    for x in Data: 
        getPeaks(wfs=x, fft_template_conj=fft_template_fit_conj,Template=Template)
    raw_input('')


def v(t, curve, time):
    #The original data, shifted down to the zxis
    # para = fit_func(time, curve)
    # set = FuncExpGausMulti(time, *para)
    baseline = 0
    data = curve - baseline
    v = np.fft.rfft(data[2750:7250])
    return v[t]

def s_star(t, curve, time, Template):
    s_star = np.conj(np.fft.rfft(Template))
    return s_star[t]

def s(t, curve, time, Template):
    s = np.fft.rfft(Template)
    return s[t]

def vary_pulse_time(curve , time, Template):
    array = []
    J = abs(np.fft.rfft(curve[np.where(time<-100)])) ** 2
    for t0_hat in np.linspace(time[0], time[-1], len(time)/5000):
        print(t0_hat)
        n = range(0, 2251, 1)
        j = complex(0, 1)
        n_pi = np.linspace(-2251 * np.pi, 2251 * np.pi, 2251)
        omega = n_pi / 2251
        shift_term = np.exp(j * t0_hat * omega)
        top = np.sum((shift_term) * s_star(n, curve, time, Template) * v(n, curve, time) / J).real
        bottom = np.sum(np.abs(s(n, curve, time, Template)) ** 2 / J)
        A = top / ( bottom)
        array = np.append(array, A)
    return array


def plot_fft(fft, zoom=False):
    for i in range(len(fft)):
        plt.plot(range(fft[i].size), fft[i])
    plt.xlabel('idx')
    plt.ylabel('fft')
    plt.grid(True)
    plt.xlim(xmin=0.0, xmax=15000)
    plt.yscale('log')
    

def plot_wfs(wfs, label, lim=[2200,4700]):
    for i in range(len(wfs)):
        plt.plot(range(len(wfs[i])), wfs[i], label=label[i])
    plt.xlabel('idx')
    plt.ylabel('voltage')
    plt.grid(True)
    plt.legend(loc="best")

def getPeaks(wfs, fft_template_conj,Template):
    wfs = np.asarray(wfs)
    wfs_empty = np.arange(wfs.size)
    peaks, peak_time, j = [], 0, 0
    fft_i = np.fft.rfft(wfs)
    conv_i = fft_i * fft_template_conj
    ifft_i = np.fft.irfft(conv_i)
    ifft_i = (ifft_i.real - np.min(ifft_i.real)) / (np.max(ifft_i.real) - np.min(ifft_i.real))

    ifft_i_base = ifft_i[abs(ifft_i) < (abs(ifft_i.mean()) + 1. * ifft_i.std())]
    thresh = 5. * ifft_i_base.std() + ifft_i_base.mean()
    # plot_wfs(wfs=[ifft_i, [thresh] * len(ifft_i), wfs, Template], label=['i-ifft_{fit}', 'threshold_{fit}', 'original', 'template'])

    plt.ion()
    plt.plot(np.arange(ifft_i.size), ifft_i)
    plt.plot(range(len(ifft_i)),[thresh] * len(ifft_i))
    plt.draw()
    

    return peaks
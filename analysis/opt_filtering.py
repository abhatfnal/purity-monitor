import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.special as spl
import glob

def main(J, FileName):
    load_data(FileName)
    cathode = np.load('cathode_' + str(FileName) + '.npy')
    anode = np.load('anode_' + str(FileName) + '.npy')
    time = np.load('time_' + str(FileName) + '.npy')
    cathode_para = fit_func(time, cathode)
    anode_para = fit_func(time, anode)
    print('Cathode Rise Time: ' + str(cathode_para[1]))
    print('Anode Rise Time: ' + str(anode_para[1]))
    cathode_fit = fit_func_normal_cathode(time, cathode)
    anode_fit = fit_func_normal_anode(time, anode)

    cathode_amp_array = vary_pulse_time(cathode, time)
    anode_amp_array = vary_pulse_time(anode, time)
    cathode_amp = np.min(cathode_amp_array)
    anode_amp = np.max(anode_amp_array)
    print cathode_amp_array
    # plt.plot(cathode)
    plt.plot(cathode_amp_array)
    plt.show()
    print('Cathode: ' + str(cathode_amp))
    print('Anode:' + str(anode_amp))
    
def fit_func(time, curve):
    #This gives a non-normalized fit of the pulse based on
    p0 = np.mean(curve[0:1000])
    p1 = 2 * 10 ** (-7)
    p2 = .000257
    p3 = 1
    p4 = -.1435
    #if curve[5] > curve[5000]:
        #p4 = -.1435
    #else:
        #p4 = .1435
    p5 = 0
    #Estimated values for p
    p = p0, p1, p2, p3, p4, p5
    yvals = curve
    xvals = time
    para, error = opt.curve_fit(FuncExpGausMulti, xdata=xvals, ydata=yvals, p0=p)
    return para

def fit_func_normal_cathode(time, curve):
    para = fit_func(time, curve)
    set = FuncExpGausMulti(time, *para)
    baseline = set[5]
    amp = np.amin(set) - baseline
    fit_func_normal = (set - baseline)/amp
    return fit_func_normal

def fit_func_normal_anode(time, curve):
    para = fit_func(time, curve)
    set = FuncExpGausMulti(time, *para)
    baseline = set[5]
    amp = np.amax(set) - baseline
    fit_func_normal = (set - baseline)/amp
    return fit_func_normal

def fit_func_normal(time, curve):
    #Normalizes the fitted curve
    para = fit_func(time, curve)
    set = FuncExpGausMulti(time, *para)
    baseline = set[5]
    if set[5] > set[5000]:
        amp = np.amin(set) - baseline
    else:
        amp = np.amax(set) - baseline
    fit_func_normal = (set - baseline)/amp
    return fit_func_normal

def load_data(FileName):
    array = open(str(FileName) + '.txt', 'r')
    new = array.read()
    new2 = np.array(new.splitlines())

    time = []
    cathode = []
    anode = []

    def get_value(string):
        exp = int(string[-2:])
        root = float(string[:-3])
        value = root * (10**exp)
        return value

    for i in range(0, 10000, 1):
        line = str(new2[i])
        data_line = line.split('\t')
        value1 = data_line[0]
        time_val = get_value(value1)
        time = np.append(time, time_val)
        value2 = data_line[1]
        cath_val = get_value(value2)
        cathode = np.append(cathode, cath_val)
        value3 = data_line[2]
        anode_val = get_value(value3)
        anode = np.append(anode, anode_val)

    np.save('cathode_' + str(FileName), cathode)
    np.save('anode_' + str(FileName), anode)
    np.save('time_' + str(FileName), time)

def FuncExpGausMulti(x, p0, p1, p2, p3, p4, p5):
#This produces the function 3.5
# // p0: baseline, B
# // p1: gaussian sig
# // p2: exponential decay constant
# // p3: number of pulses
# // p[4+2*i]: pulse[i] amplitude
    x = np.asarray(x)
    val = np.full((x.size,),p0)
    time = x-p5
    #This is the term time = (t - mu)/tau
    val_tot = val + (p4/2.* np.exp((p1*p1/p2/2.-time)/p2)* spl.erfc((p1*p1/p2-time)/np.sqrt(2)/p1))
            #This is the function is defined in eqn 3.5                                                                        )
    return val_tot

def v(t, curve, time):
    #The original data, shifted down to the zxis
    para = fit_func(time, curve)
    set = FuncExpGausMulti(time, *para)
    baseline = set[5]
    data = curve - baseline
    v = np.fft.rfft(data[2750:7250])
    return v[t]

def s_star(t, curve, time):
    s_star = np.conj(np.fft.rfft(fit_func_normal(time, curve)[2750:7250]))
    return s_star[t]

def s(t, curve, time):
    s = np.fft.rfft(fit_func_normal(time, curve))
    return s[t]

def vary_pulse_time(curve , time):
    array = []
    for t0_hat in np.linspace(time[0], time[-1], len(time)/2000):
        n = range(0, 2251, 1)
        j = complex(0, 1)
        n_pi = np.linspace(-2251 * np.pi, 2251 * np.pi, 2251)
        omega = n_pi / 2251
        shift_term = np.exp(j * t0_hat * omega)
        top = np.sum((shift_term) * s_star(n, curve, time) * v(n, curve, time) / J).real
        bottom = np.sum(np.abs(s(n, curve, time)) ** 2 / J)
        A = top / ( bottom)
        array = np.append(array, A)
    return array


J = np.load('varying_volt_j_anode.npy')
main(J, 'field_data_1')
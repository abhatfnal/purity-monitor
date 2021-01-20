import matplotlib
import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

startTime = datetime.datetime.now()

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex' : True,'font.size' : 14,'font.family' : 'lmodern','text.latex.unicode': True,}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']

def GetData(Data, Filename, Cut=0): 
    Data['Date'] = []
    Data['Temperature 1'] = []
    Data['Temperature 2'] = []
    Data['Humidity'] = []

    if Cut == -1: 
        Cut = 1E10

    File = open(Filename, 'r')
    for ii,line in enumerate(File):
        values = line.split('\t')
        dateval = datetime.datetime.strptime(values[1]+' '+values[2], '%b %d %Y %H:%M:%S')
        Data['Date'].append(dateval)
        Data['Temperature 1'].append(float(values[3]))
        Data['Humidity'].append(float(values[4]))
        Data['Temperature 2'].append(float(values[5]))
        if ii % 50000 == 0: 
            print(ii)
        if ii > Cut:
            break
    
    return Data

def PlotData(Data): 
    fig, ax1 = plt.subplots(figsize=(12,6))

    plt.gcf().autofmt_xdate()
    formatter = matplotlib.dates.DateFormatter('%b-%d %H:%M')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    ax1.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=12))
    ax1.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(6))
    ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    ax1.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))
    ax1.grid(b=True, which='major', color='k', linestyle='--', alpha=0.8)
    ax1.grid(b=True, which='minor', color='grey', linestyle=':')

    plt.xlabel('Time[mm/dd  HH:MM]', fontsize=14)
    plt.ylabel(r'Temperature [$^\circ$C]', fontsize=14)
    plt.xlim(np.min(Data['Date']), np.max(Data['Date']))
    plt.ylim(23, 25)

    TemperatureLabel1 = r'Temperature 1: $%.2f \pm %.2f$' % (np.mean(Data['Temperature 1']), np.std(Data['Temperature 1']))
    TemperatureLabel2 = r'Temperature 2: $%.2f \pm %.2f$' % (np.mean(Data['Temperature 2']), np.std(Data['Temperature 2']))
    HumidityLabel = r'Humidity: $%.2f \pm %.2f$' % (np.mean(Data['Humidity']), np.std(Data['Humidity']))

    ax1.plot(Data['Date'], Data['Temperature 1'], label=TemperatureLabel1, color = colors[0], linewidth=1)
    ax1.plot(Data['Date'], Data['Temperature 2'], label=TemperatureLabel2, color = colors[1], linewidth=1)
    plt.legend(loc='upper left')

    ax2 = ax1.twinx()
    formatter = matplotlib.dates.DateFormatter('%m/%d  %H:%M')
    plt.gcf().axes[1].xaxis.set_major_formatter(formatter)

    ax2.plot(Data['Date'], Data['Humidity'], label=HumidityLabel, color = colors[2], linewidth=1)
    plt.ylim(20, 45)
    plt.ylabel(r'Humidity [$\%$]', fontsize=14)
    plt.legend(loc='upper right')

    fig.tight_layout()
    plt.savefig('sensor_data.pdf')


def PlotFFT(Data):
    import numpy.fft as fft
    spectrum = fft.fft(Data['Temperature 1'])
    freq = fft.fftfreq(len(spectrum))

    fig, ax1 = plt.subplots(figsize=(12,6))
    plt.xlabel('Period [Hours]', fontsize=14)
    plt.ylabel('Amplitude', fontsize=14)
    plt.plot((1.0/abs(freq))/3600.0, abs(spectrum))
    plt.grid()

    plt.xscale('log')
    plt.yscale('log')
    fig.tight_layout()
    plt.savefig('sensor_data_fft.pdf')
    


if __name__ == '__main__':
    Filename = '/project/fas/david_moore/aj487/purity_monitor/data2.txt'
    Data = {}
    Data = GetData(Data = Data, Filename = Filename, Cut = -1)
    PlotData(Data)
    PlotFFT(Data)
    # plt.show()
    plt.close()
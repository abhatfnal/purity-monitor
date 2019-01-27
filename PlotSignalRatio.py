import matplotlib.pyplot as plt
from uncertainties import ufloat, nominal_value, std_dev
from matplotlib.ticker import MultipleLocator

filename = 'data/LiquidXenonSignalRatio.dat'
filename2 = 'data/LiquidXenonSignalRatio2.dat'

def Run(filename):
    drift = []
    anode = []
    cathode = []
    anodeMax = []
    cathodeMax = []

    file = open(filename, "r")
    for ii, line in enumerate(file):
        if(ii==0):
            continue
        columns = line.split()
        drift.append(float(columns[0]))
        anode.append(ufloat(columns[1], columns[3]))
        cathode.append(ufloat(columns[2], columns[4]))
        anodeMax.append(float(columns[5]))
        cathodeMax.append(float(columns[6]))

    ratio = [nominal_value(x/y) for x,y in zip(anode,cathode)]
    ratioErr = [std_dev(x/y) for x,y in zip(anode,cathode)]
    ratioMax = [x/y for x,y in zip(anodeMax,cathodeMax)]

    fig = plt.figure(figsize=(12,9))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Drift Length [cm]', fontsize=14)
    plt.ylabel('Anode/Cathode', fontsize=14)

    plt.xlim(0, 16)
    plt.ylim(0.0, 1.2)
    plt.errorbar(drift, ratio, yerr=ratioErr, fmt='X', label='Fit')
    plt.scatter(drift, ratioMax, color='r', label='Max')
    plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)
    plt.savefig("plots/LiquidXenonSignalRatio.pdf", bbox_inches = 'tight', pad_inches = 0.2)
    plt.show()


for f in [filename]:
    Run(f)

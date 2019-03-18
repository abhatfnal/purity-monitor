import matplotlib.pyplot as plt
from uncertainties import ufloat, nominal_value, std_dev
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
import numpy as np


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


filename = 'data/20190131_drift.dat'
file = open(filename, "r")

amp = []
drift = []
dist = []



for ii, line in enumerate(file):
    if(ii==0):
        continue
    columns = line.split()
    dist.append(float(columns[0]))
    drift.append(float(columns[1]))
    amp.append(float(columns[2]))

y = func(np.asarray(dist), 100, 1, 1)
popt, pcov = curve_fit(func, dist, amp)
fig = plt.figure(figsize=(12,9))
ax = fig.gca()
ax.grid(b=True, which='major', color='k', linestyle='--')
ax.grid(b=True, which='minor', color='grey', linestyle=':')
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(5))

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlabel('Drift Length [cm]', fontsize=14)
plt.ylabel('Anode/Cathode', fontsize=14)

# plt.xlim(0, 16)
plt.ylim(0.0, 200)
plt.scatter(dist, amp, color='r', label='Data')
plt.plot(dist, func(dist, *popt), 'r-', label='fit')
plt.legend(loc='upper right', prop={'size': 14}, numpoints=1)
plt.savefig("plots/LiquidXenonSignalRatio.pdf", bbox_inches = 'tight', pad_inches = 0.2)
plt.show()

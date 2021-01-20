import numpy as np
import matplotlib
import os, time, datetime
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator
from decimal import Decimal
from scipy.optimize import curve_fit

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']
parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath", nargs="*")
arg = parser.parse_args() 

def func(x,amp,alpha,epsilon):
    return amp * alpha * x * ((x + (1/alpha) ) / (x + epsilon) )
    # return amp * x * ((x -55.01 ) / (x + epsilon) )

if __name__ == '__main__':
    fig = plt.figure(figsize=(12,9))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Flow Rate [SCCPS]', fontsize=14)
    plt.ylabel('Life Time [us]', fontsize=14)

    for jj, x in enumerate(arg.filepath): 
        flow = []
        lifetime = []
        file = open(x, "r")
        path, filename = os.path.split(x)
        for ii, line in enumerate(file):
            flow.append(float(line[0:6]))
            lifetime.append(float(line[8:]))
            print float(line[0:6]), float(line[8:])
        flow = np.array(flow)
        lifetime = np.array(lifetime)
        plt.scatter(flow, lifetime, label=filename, linewidth=1.5, c=colors[0], edgecolor='face')

        x = np.linspace(0,350,350)
        init_vals = [1,1,1]
        best_vals, covar = curve_fit(func, flow, lifetime, p0=init_vals)
        print best_vals

        plt.plot(x, func(x, best_vals[0], best_vals[1], best_vals[2]), label=r'f(x)=$%f \cdot %d \cdot x \cdot \frac{x+1/%d}{x+%d}$' % (best_vals[0], best_vals[1], best_vals[1], best_vals[2]), linewidth=1.5, c=colors[1])
        ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0, fontsize=12)
        plt.xlim(0,365)  
        plt.ylim(0,11)
        plt.show()

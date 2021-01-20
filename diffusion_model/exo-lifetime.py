import numpy as np
import matplotlib
import os, time, datetime
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import MultipleLocator
from decimal import Decimal
from scipy.optimize import curve_fit

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 12,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']
parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, action="store", dest="filepath", nargs="*")
arg = parser.parse_args() 

def func(x,amp,alpha,epsilon):
    return amp * alpha * x * ((x + (1/alpha) ) / (x + epsilon) )
    # return amp * x * ((x -55.01 ) / (x + epsilon) )

if __name__ == '__main__':
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    # ax.xaxis.set_major_locator(MultipleLocator(1000))
    # ax.xaxis.set_minor_locator(MultipleLocator(200))
    # ax.yaxis.set_major_locator(MultipleLocator(1000))
    # ax.yaxis.set_minor_locator(MultipleLocator(200))

    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlabel('Run Number', fontsize=14)
    plt.ylabel(r'Life Time [$\mu$s]', fontsize=14)

    for jj, x in enumerate(arg.filepath): 
        run = []
        lifetime = []
        errorlow = []
        errorhigh = []
        file = open(x, "r")
        path, filename = os.path.split(x)
        for ii, line in enumerate(file):
            data = line.split(" ")
            if(float(data[6])>0 and float(data[6])<10000):
                run.append(float(data[1]))
                lifetime.append(float(data[6]))
                errorlow.append(float(data[7]))
                errorhigh.append(float(data[8]))
        
        errorlow = np.array(errorlow)
        errorhigh = np.array(errorhigh)
        run = np.array(run)
        lifetime = np.array(lifetime)

        print np.mean(lifetime[np.where(run<6700)[0]]), np.std(lifetime[np.where(run<6700)[0]])
        print np.mean(lifetime[np.where(run>6700)[0]]), np.std(lifetime[np.where(run>6700)[0]])

        x = np.linspace(0,np.max(run),np.max(run))
        y = [np.mean(lifetime[np.where(run<6700)[0]])]*len(x)
        y2 = [np.mean(lifetime[np.where(run>6700)[0]])]*len(x)
        combined = np.vstack((errorlow, errorhigh))
        plt.scatter(run, lifetime, label='', linewidth=1.0, s=4.0, c=colors[0], edgecolor='face')
        # plt.errorbar(run, lifetime, yerr=combined, linestyle="None", c=colors[0])
        # x = np.linspace(0,350,350)
        # init_vals = [1,1,1]
        # best_vals, covar = curve_fit(func, flow, lifetime, p0=init_vals)
        # print best_vals
        # plt.plot(x, func(x, best_vals[0], best_vals[1], best_vals[2]), label=r'f(x)=$%f \cdot %d \cdot x \cdot \frac{x+1/%d}{x+%d}$' % (best_vals[0], best_vals[1], best_vals[1], best_vals[2]), linewidth=1.5, c=colors[1])

        
        ax.legend(loc='upper right', fontsize=12)
        # plt.xlim(0,365)  
        plt.ylim(0,7000) 
        # fig.tight_layout()
        plt.savefig('exo-lifetime.pdf', dpi=1000, bbox_inches='tight')
        plt.show()

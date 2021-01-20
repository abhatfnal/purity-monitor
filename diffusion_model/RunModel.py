import numpy as np
import matplotlib.pyplot as plt
import Library as Lib 
import Outgassing as Out
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, LogLocator

params = {'text.usetex':True,'font.size':12,'font.family':'serif'}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']


class System(): 
    def __init__(self, SystemName, Material, Solute, Dimension):
        self.Name = SystemName
        self.Material = Material
        self.Solute = Solute
        self.Dimension = Dimension
        self.XeMass = Lib.System.get(SystemName).get('Xenon Mass')
        
        self.Diffusion = Lib.Material.get(Material).get(Solute).get('Diffusion Constant')
        self.Solubility = Lib.Material.get(Material).get(Solute).get('Solubility')
        self.ActivationEnergy = Lib.Material.get(Material).get(Solute).get('Activation Energy')
        self.Abundance =  Lib.Gas.get(Solute).get('Abundance in Air')
        self.MolarMass =  Lib.Gas.get(Solute).get('Molar Mass')

        self.Volume = Lib.System.get(SystemName).get(Material).get(Dimension).get('Volume')
        self.Area = Lib.System.get(SystemName).get(Material).get(Dimension).get('Area')
        self.Thickness = Lib.System.get(SystemName).get(Material).get(Dimension).get('Thickness')

    def Print(self):
        Attributes = vars(self)
        print(', '.join("%s: %s" % item for item in Attributes.items()))


def GetTimeStamps(Points, Spacing, TimeScale): 
    X = []
    for ii, x in enumerate(Points): 
        if ii == len(Points)-1: break 
        X.append(np.linspace(Points[ii], Points[ii+1], int((Points[ii+1] - Points[ii]) / Spacing + 1)))
    return np.array(X)

def PlotImpuritiesVsTime(Data):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    plt.xlabel('Time [hours]', fontsize=16)
    plt.ylabel(r'Total Number of Impurities', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(LogLocator(base=10,numticks=12,))
    for jj,data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.Impurities)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[jj], linewidth=2.)

    plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    # plt.ylim(ymin=1E12,ymax=1E20)
    ax.legend(loc='upper right', fontsize=14)
    fig.tight_layout()

def PlotFlowRateVsTime(Data):
    fig = plt.figure(figsize=(10,7))
    ax = fig.gca()
    plt.xlabel('Time [Hours]', fontsize=16)
    plt.ylabel(r'Outgassing Rate [mBar$\,\cdot\,$Liter/s]', fontsize=16)
    plt.yscale('log')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    ax.grid(b=True, which='major', color='grey', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(LogLocator(base=10,numticks=12,))
    for jj,data in enumerate(Data):
        for ii,(X,Y) in enumerate(zip(data.Time, data.FlowRate)):
            plt.plot(X, Y, label=data.Labels[ii], color=colors[jj], linewidth=2.)

    plt.xlim(np.min(Data[0].Time[0]), np.max(Data[0].Time[-1]))
    plt.ylim(ymin=1E-12)
    ax.legend(loc='upper right', fontsize=14)
    fig.tight_layout()

if __name__ == '__main__':

    EXO = System('YLXPS', 'Teflon', 'Oxygen', 'EXO-Teflon')
    EXO2 = System('YLXPS', 'Teflon', 'Oxygen', 'Stock-Teflon')
    EXO3 = System('YLXPS', 'Teflon', 'Oxygen', 'Columbia-Teflon')

    Systems = [EXO,EXO2,EXO3]
    Labels = [['EXO-200 Teflon'], ['Stock Room Teflon'], ['Columbia Setup Teflon']]

    for ii,System in enumerate(Systems): 
        System.Temp = [295,295]       
        System.InitialImpurities = Out.GetInitialImpurities(System, '#')
        System.DiffConstants = Out.GetDiffTemp(System, Temperatures=System.Temp)
        System.Labels = Labels[ii]
        System.Time = GetTimeStamps(Points=[0,100], Spacing=0.1, TimeScale='Hours')
        System.Impurities = Out.GetImpuritiesVsTime(Data=System, TimeScale='Hours')
        System.FlowRate = Out.GetFlowRateVsTime(Data=System, Units='mBar Liter', TimeScale='Hours')
    PlotImpuritiesVsTime(Systems)
    plt.savefig('impurities.pdf')
    PlotFlowRateVsTime(Systems)
    plt.savefig('outgassing_rate.pdf')
    plt.show()

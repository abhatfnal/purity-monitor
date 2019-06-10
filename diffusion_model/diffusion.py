import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.misc import factorial as fact
from decimal import Decimal
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 12,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']


path = '/home/fas/david_moore/aj487/purity_monitor/plots/diffusion_model/'

def GetConcentration(x, Diff, Thickness, Conc): 
    value = 0.0 
    for N in range(0,1000,1): 
        factor1 = 1.0/((2.0*N+1.0)**2)
        factor2 = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        comp = factor1*factor2
        value = value + comp
    value = value * (Conc*8.0*Thickness)/(np.pi**2*2.0)
    return value 

def GetFlowRate(x, Diff, Thickness, Conc, Area): 
    value = 0.0 
    Conc = Conc*Area
    for N in range(0,1000,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        value = value + factor
    value = value * (Conc*8.0)/Thickness 
    return value

def GetDiffTemp(Diff, Temp, Energy=0.17): 
    kb = 8.6173303E-5 # eV/K
    DiffTemp = Diff * np.exp(Energy/kb * ((1.0/294.15) - (1.0/Temp)))
    print ' | Diffusion coefficient at %s K:\t' % Temp, DiffTemp 
    return DiffTemp

def PlotDiffVsTemp():
    ActEn = [] 
    for ii,x in enumerate(range(20)):
        ActEn.append(0.05*(ii+1))
    fig = plt.figure(figsize=(12,7))
    ax = fig.gca()
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(5))

    cm = plt.get_cmap('gist_rainbow')
    NUM_COLORS = 20
    ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.xlabel('Temperature [K]', fontsize=16)
    plt.ylabel(r'Diffusion Coefficient [$\mathrm{cm}^2 \mathrm{s}^{-1}$]', fontsize=16)
    # plt.xscale('log')
    plt.yscale('log')
    x = np.linspace(150,400,400-150)
    for val in ActEn:
        Diff0 = Diff/np.exp(-val/(kb*(21+273.15)))
        y = Diff0/timetodays*np.exp(-val/(kb*x))
        plt.plot(x, y, label='%.2f eV'%val, linewidth=1.5)
    Diff0 = Diff/np.exp(-0.17/(kb*(21+273.15)))
    y = Diff0/timetodays*np.exp(-0.17/(kb*x))
    plt.plot(x, y, label='0.17 eV', color='black', linewidth=1.5, linestyle='--')
    ax.legend(loc='lower left', bbox_to_anchor=(1, 0), ncol=1, borderaxespad=0, fontsize=10)
    
def GetActEnergy(StartEnergy, NumEnergy, StepEnergy):
    ActEn = [] 
    for ii in range(NumEnergy):
        ActEn.append(StartEnergy + StepEnergy*(ii+1))
    return ActEn

def GetLegend(LegendData, Units):
    Legend = []
    for x in LegendData:
        Legend.append((r'%.2f$\,$'+Units)% x)
    return Legend

def PlotEXO200Teflon():
    thickness = 0.15 # cm
    timetodays = 3600*24
    Ea = 0.17 #eV
    C0 = 3E5 # ppt exo-200 
    Diff = 31.4E-8*timetodays # cm^2/day

    days = 0.1
    days2 = days+1
    days3 = days2+3000
    multi = 100
    xall = []
    x = np.linspace(0.0, days, days*multi+1)
    x2 = np.linspace(days, days2, (days2-days)*multi+1)
    x3 = np.linspace(days2, days3, (days3-days2)*multi+1)
    xs = np.append(x,x2)
    xs = np.append(xs,x3)
    xall.append(x)
    xall.append(x2)
    xall.append(x3)

    ActEn = GetActEnergy(0.09, 15, 0.03)
    ActEn = [0.17]
    yall = []
    ys = [] 
    xss = []
    for ii,val in enumerate(ActEn):
        y = GetConcentration(x, Diff, thickness, C0)
        yall.append(y)
        y2 = GetConcentration(x2, GetDiffTemp(Diff, Temp=273.15+50.0, Energy=val), thickness, y[-1])
        y2 = y2 * y[-1]/y2[0]
        y = np.append(y, y2)
        yall.append(y2)
        y3 = GetConcentration(x3, GetDiffTemp(Diff, Temp=161.0,Energy=val), thickness, y2[-1])
        y3 = y3 * y2[-1]/y3[0] 
        y = np.append(y, y3)
        ys.append(y)
        yall.append(y3)
        xss.append(xs)

    # Legend = GetLegend(ActEn, 'eV')
    Legend = ['Pump Down', 'Warm Xenon Gas', 'Liquid Xenon']
    XRange = [1E-2,np.max(xs)]
    YRange = [1E-30,np.max(ys)*10]
    PlotVsTime(XValue=xall,YValue=yall,Legend=Legend,XRange=XRange,YRange=YRange)
    plt.savefig(path+'exo_teflon_conc.pdf', dpi=1000, bbox_inches='tight')
    # plt.show()

    Label = ['Time [days]', r'Gas Flow Rate [$\mathrm{Molecules}\,\mathrm{day}^{-1}$]']
    C00 = 0.046/(693*32) * 6.023E23 
    Area = 2 * np.pi * 18.3 * 40.0
    C00 = C00 * Area
    ys = [] 
    yall = []
    for ii,val in enumerate(ActEn):
        y = GetFlowRate(x, Diff, thickness, C0, Area)
        yall.append(y)
        y2 = GetFlowRate(x2, GetDiffTemp(Diff, 273.15+50,val), thickness, y[-1], Area)
        y2 = y2 * y[-1]/y2[0]
        y = np.append(y, y2)
        yall.append(y2)
        y3 = GetFlowRate(x3, GetDiffTemp(Diff, 161.0,val), thickness, y2[-1], Area)
        y3 = y3 * y2[-1]/y3[0]
        y = np.append(y, y3)
        ys.append(y)
        yall.append(y3)
    XRange = [1E-2,np.max(xs)]
    YRange = [1E-30,np.max(ys)*10]
    PlotVsTime(XValue=xall,YValue=yall,Legend=Legend,Label=Label,XRange=XRange,YRange=YRange)
    plt.savefig(path+'exo_teflon_flow.pdf', dpi=1000, bbox_inches='tight')
    # plt.show()

def PlotnEXOKapton():
    C0Kapton = 11156 # ppt nEXO 
    DKapton = 31.4E-8 # cm^2/day
    TKapton = 0.005 # cm
    days = 100000
    x = np.linspace(0.0, days, days*1+1)

    y = GetConcentration(x, DKapton, TKapton, C0Kapton)
    y2 = GetConcentration(x, GetDiffTemp(Diff=DKapton, Temp=161.0), TKapton, C0Kapton)
    Legend = ['20 C','-110 C']
    Label = ['Time [sec]','Impurity Concentration [ppt]']
    XRange = [1E-1,np.max(x)]
    PlotVsTime(XValue=[x,x,x],YValue=[y,y2],Legend=Legend,XRange=XRange,Label=Label)
    plt.savefig(path+'nexo_kapton_conc.pdf', dpi=1000, bbox_inches='tight')

    Area = 12500.0 #cm^2
    C0KaptonN = 1.531E+14
    y = GetFlowRate(x, DKapton, TKapton, C0KaptonN, Area)
    y2 = GetFlowRate(x, GetDiffTemp(DKapton, Temp=161.0), TKapton, C0KaptonN, Area)
    Label = ['Time [sec]',r'Gas Flow Rate [$\mathrm{Molecules}\,\mathrm{sec}^{-1}$]']
    PlotVsTime(XValue=[x,x],YValue=[y,y2],Legend=Legend,XRange=XRange,Label=Label)
    plt.savefig(path+'nexo_kapton_flow.pdf', dpi=1000, bbox_inches='tight')

    # plt.show()

def PlotVsTime(XValue,YValue,Legend,Title='',XRange=0,YRange=0,Ticks=[50,10],Label=['',''],Log=[True,True],LPos=(0,0),Size=(6,5)):
    fig = plt.figure(figsize=Size)
    ax = fig.gca()
    # color = ['blue','red','green','black','purple','orange']
    Xtot = ()
    Ytot = ()
    for x,y in zip(XValue,YValue):
        Xtot = np.append(Xtot,x)
        Ytot = np.append(Ytot,y)
    if(Label==['','']):
        plt.xlabel('Time [days]', fontsize=16)
        plt.ylabel(r'Impurity Concentration [ppt]', fontsize=16)
    else:
        plt.xlabel(Label[0], fontsize=16)
        plt.ylabel(Label[1], fontsize=16)
    if(Log[0]): 
        plt.xscale('log')
    if(Log[1]): 
        plt.yscale('log')
    if(XRange==0):
        plt.xlim(0.1,np.max(Xtot))
    else: 
        plt.xlim(XRange[0],XRange[1])
    if(YRange==0):
        plt.ylim(1E-10,np.max(Ytot)*10)
    else: 
        plt.ylim(YRange[0],YRange[1])
    for ii, (xval,yval) in enumerate(zip(XValue,YValue)):
        plt.plot(xval, yval, color=colors[ii], label=Legend[ii], linewidth=1.5)
    if(LPos!=(0,0)):
        ax.legend(loc='lower left', bbox_to_anchor=LPos, ncol=1, borderaxespad=0, fontsize=14)
    else:
        ax.legend(loc='upper right',fontsize=16)
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    # ax.xaxis.set_major_locator(MultipleLocator(Ticks[0]))
    # ax.xaxis.set_minor_locator(MultipleLocator(Ticks[1]))
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.title(Title)
    fig.tight_layout()

def PlotLZTeflon():
    timetodays = 3600*24
    thickness = 2.0
    #oxygen
    C0 = 943000.0 # ppt
    Diff = 31.4E-8*timetodays # cm^2/day
    #nitrogen
    C02 = 1491000.0 # ppt
    Diff2 = 15.1E-8*timetodays # cm^2/day
    #krypton
    C03 = 290.5 # ppt
    Diff3 = 5.6E-8*timetodays # cm^2/day

    days = 500
    multi = 100
    x = np.linspace(0.0, days, days*multi+1)
    y = GetConcentration(x, Diff, thickness, C0)
    y2 = GetConcentration(x, Diff2, thickness, C02)
    y3 = GetConcentration(x, Diff3, thickness, C03)

    YRange = [1E-2,1E7]
    Label = ['Time [days]', r'Impurity Concentration [ppt]']
    Legend = ['Oxygen', 'Nitrogen', 'Krypton']
    Log = [False,True]
    PlotVsTime(XValue=[x,x,x],YValue=[y,y2,y3],Legend=Legend,Label=Label,Log=Log,YRange=YRange,Size=(10,5))
    # plt.show()
    plt.savefig(path+'lz_teflon_conc.pdf', dpi=1000, bbox_inches='tight')
    # plt.savefig(path+'example.pgf')

def PlotEXOSteel():
    sec = 5000000
    multi = 1
    Flow = 2E-7 *60 #mbar*liter/s
    Thickness = 0.165 #cm 
    DOxygen = 31.4E-8
    Diff = GetDiffTemp(DOxygen, Temp=161.0)
    x = np.linspace(0.0, sec, sec*multi+1)
    y = 0
    y2 = 0 
    for N in range(0,10,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        y = y + factor
    y = y * Flow
    y2 = 0 
    for N in range(0,50,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        y2 = y2 + factor
    y2 = y2 * Flow
    y3 = 0 
    for N in range(0,100,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        y3 = y3 + factor
    y3 = y3 * Flow
    Log = [True,True]
    Legend = ['10','50','100']
    Label = ['Time [sec]',r'Gas Flow Rate [$\mathrm{mbar}\,\mathrm{liter}\,\mathrm{sec}^{-1}$]']
    XRange = [1E-1,np.max(x)]
    YRange = [1E-7,1E0]

    PlotVsTime(XValue=[x,x,x],YValue=[y,y2,y3],Legend=Legend,Label=Label,Log=Log,XRange=XRange,YRange=YRange,Size=(10,10))
    plt.show()
    print np.sum(y[0:60*24])

if __name__ == '__main__':
    # PlotDiffVsTemp()
    # PlotnEXOKapton()
    # PlotEXO200Teflon()
    # PlotLZTeflon()
    PlotEXOSteel()

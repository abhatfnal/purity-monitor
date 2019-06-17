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

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']


path = '/home/fas/david_moore/aj487/purity_monitor/plots/diffusion_model/'
D0 = 0.000225323


def GetConcentration(x, Diff, Thickness, Conc): 
    value = 0.0 
    for N in range(0,2,1): 
        factor1 = 1.0/((2.0*N+1.0)**2)
        factor2 = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        comp = factor1*factor2 * (Conc*8.0*Thickness)/(np.pi**2*2.0)
        value = value + comp
    return value 

def GetFlowRate(x, Diff, Thickness, Conc, Area): 
    value = 0.0 
    for N in range(0,2,1): 
        factor = np.exp(-(np.pi*(2.0*N+1.0)/Thickness)**2*Diff*x)
        value = value + factor*(8.0*Conc*Diff)/Thickness
    value = value * Area 
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

    # cm = plt.get_cmap('gist_rainbow')
    # NUM_COLORS = 20
    # ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    colors = ['#415E80', '#3694FF', '82BCFF', '003A7D', '6896CB', 'B5D6FF', '', '', '', '', '', '', '', '', '', '', '', '', '',]
    palette = plt.get_cmap('Set1')

    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.xlabel('Temperature [K]', fontsize=16)
    plt.ylabel(r'Diffusion Coefficient [$\mathrm{cm}^2 \mathrm{s}^{-1}$]', fontsize=16)
    # plt.xscale('log')
    plt.yscale('log')
    x = np.linspace(150,400,400-150)
    for ii,val in enumerate(ActEn):
        Diff0 = Diff/np.exp(-val/(kb*(21+273.15)))
        y = Diff0/timetodays*np.exp(-val/(kb*x))
        plt.plot(x, y, label='%.2f eV'%val, linewidth=1.5, color=palette(num))
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

def GetXValues(days, multi):
    xs = []
    xall = []
    xind = []
    for ii,day in enumerate(days): 
        if ii==0:
            x = np.linspace(0.0, days[ii], days[ii]*multi[ii]+1)
            xind.append(np.linspace(0.0, days[ii], days[ii]*multi[ii]+1))
        else:
            x = np.linspace(days[ii-1], days[ii], (days[ii]-days[ii-1])*multi[ii]+1)
            xind.append(np.linspace(0.0, days[ii]-days[ii-1], (days[ii]-days[ii-1])*multi[ii]+1))
        xs = np.append(xs,x)
        xall.append(x)
    return np.asarray([xs]), xall, xind

def PlotEXO200Teflon():
    thickness = 0.15 # cm
    Ea = 0.17 #eV
    C0 = 8.6E20 # number of total impurities in teflon 
    Diffs = [0.0271296, 0.0496, 0.000105013] # cm^2/day
    days = [1,3,10]
    multi = [20, 20, 20]
    xs, xall, xind = GetXValues(days=days,multi=multi)

    PDrop = 1000.0
    PurOutput = 1E16
    y1 = GetConcentration(xind[0], Diffs[0], thickness, C0)
    for ii,yy in enumerate(y1): 
        if y1[ii]<y1[0]/PDrop: 
            y1[ii] = y1[0]/PDrop
    y2 = GetConcentration(xind[1], Diffs[1], thickness, y1[-1])
    y2 = y2 * y1[-1]/y2[0]
    for ii,yy in enumerate(y2): 
        if y2[ii]<PurOutput: 
            y2[ii] = PurOutput 
    y3 = GetConcentration(xind[2], Diffs[2], thickness, y2[-1])
    y3 = y3 * y2[-1]/y3[0]
    y4 = np.array([1.6E17]*len(xs[0]))
    y5 = np.array([9E16]*xs[0])

    index1 = np.where(y1==y1[0]/PDrop)[0][0]
    index2 = np.where(y2==PurOutput)[0][0]

    valx = [item for sublist in [xall,np.array(xs),np.array(xs)] for item in sublist]
    Legend = ['Pump Down', 'Warm Xenon Gas', 'Liquid Xenon', 'EXO-200', 'Purified']
    XRange = [1.0/np.max(multi),np.max(xs)]
    YRange = [np.min(y5[1:])/10,np.max(y1)*10]
    Label = ['Time [days]','Total Number of Impurities Left in Plastic']
    PlotVsTime(XValue=valx,YValue=[y1,y2,y3,y4,y5],XRange=XRange,YRange=YRange,Label=Label,Legend=Legend,Log=[False,True])
    plt.savefig(path+'exo_teflon_conc.pdf', dpi=1000, bbox_inches='tight')
    plt.show()

    C0 = 2.29E5*136.0/32.0       # ppt in mol/mol exo-200 
    C00 = C0 / 136           # ppt in mol/g 
    C00 = C00 * 3.1           # ppt in mol/cm^3
    C00 = C00 * 6.023E23      # ppt in #/cm^3
    C00 = C00 * 1E-12         # units of #/cm^3
    print C00

    c1 = y1/1E27*1E12 #impurity concentration ppt [mol/mol]
    c2 = y2/1E27*1E12 #impurity concentration ppt [mol/mol]
    c3 = y3/1E27*1E12 #impurity concentration ppt [mol/mol]
    c4 = y4/1E27*1E12 #impurity concentration ppt [mol/mol]
    
    YRange = [np.min(c3)/10,np.max(c1)*10]
    Label = ['Time [days]','Impurity Concentration [ppt]']
    PlotVsTime(XValue=valx,YValue=[c1,c2,c3,c4],XRange=XRange,YRange=YRange,Label=Label,Legend=Legend,Log=[False,True])
    plt.savefig(path+'exo_teflon_impur.pdf', dpi=1000, bbox_inches='tight')
    
    c1 = y1/64500 # nr of impurities/cm^3
    c2 = y2/64500 # nr of impurities/cm^3
    c3 = y3/64500 # nr of impurities/cm^3
    c4 = y4/64500 # nr of impurities/cm^3

    Area = 2 * np.pi * 18.3 * 40.0 * 2 
    f1 = GetFlowRate(xind[0], Diffs[0], thickness, c1[0], Area)
    for ii,ff in enumerate(f1): 
        if ii>=index1:
            f1[ii] = 0
    f2 = GetFlowRate(xind[1], Diffs[1], thickness, c2[0], Area)
    for ii,ff in enumerate(f2): 
        if ii>=index2:
            f2[ii] = 0
    f3 = GetFlowRate(xind[2], Diffs[2], thickness, c3[0], Area)
    f4 = [1.4E17]*len(xs[0])
    print C00, c1[-1], c2[-1]
    XRange = [1E-2,np.max(xs)]
    YRange = [np.min(f3)/10,np.max(f1)*10]
    Label = ['Time [days]', r'Outgassing Rate [$\mathrm{Molecules}\,\mathrm{day}^{-1}$]']
    PlotVsTime(XValue=valx,YValue=[f1,f2,f3,f4],Legend=Legend,Label=Label,XRange=XRange,YRange=YRange,Log=[False,True])
    plt.savefig(path+'exo_teflon_flow.pdf', dpi=1000, bbox_inches='tight')
    plt.show()

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
    plt.savefig(path+'nexo_kapton_flow.pdf', dpi=1000, bbox_inches='tight', bbox_extra_artists=(lgd))

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
        yval = [float('nan') if x==0 else x for x in yval]
        plt.plot(xval, yval, color=colors[ii], label=Legend[ii], linewidth=2.)
    if(LPos!=(0,0)):
        lgd = ax.legend(loc='upper left', bbox_to_anchor=LPos, ncol=1, borderaxespad=0, fontsize=14)
    else:
        lgd = ax.legend(loc='upper right',fontsize=16)
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='grey', linestyle=':')
    # ax.xaxis.set_major_locator(MultipleLocator(Ticks[0]))
    # ax.xaxis.set_minor_locator(MultipleLocator(Ticks[1]))
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.title(Title)
    # fig.tight_layout()

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
    plt.show()
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
    PlotEXO200Teflon()
    # PlotLZTeflon()
    # PlotEXOSteel()

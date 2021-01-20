import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from decimal import Decimal
import diffusion as dfs 

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex' : True,'font.size' : 12,'font.family' : 'lmodern','text.latex.unicode': True}
plt.rcParams.update(params) 
colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3', 'black']
path = '/home/fas/david_moore/aj487/purity_monitor/plots/diffusion_model/ylxpm/'

def PlotEXO200Teflon():
    thickness = 0.5# cm
    Ea = 0.17 #eV
    C0 = 2.6E18 # number of total impurities in teflon 
    Diffs = [0.0271296, 0.0496, 0.000105013] # cm^2/day
    Diffs = [0.0271296, 0.00000105013, 0.00000105013] # cm^2/day
    days = [4,5]
    multi = [20, 20]
    xs, xall, xind = dfs.GetXValues(days=days,multi=multi)

    PDrop = 1E10
    PurOutput = 1E16
    y1 = dfs.GetConcentration(xind[0], Diffs[0], thickness, C0)
    for ii,yy in enumerate(y1): 
        if y1[ii]<y1[0]/PDrop: 
            y1[ii] = y1[0]/PDrop
    y2 = dfs.GetConcentration(xind[1], Diffs[1], thickness, y1[-1])
    y2 = y2 * y1[-1]/y2[0]
    y4 = np.array([0.99E15]*len(xs[0]))

    index1 = 1E10
    index2 = 1E10

    valx = [item for sublist in [xall,np.array(xs),np.array(xs)] for item in sublist]
    Legend = ['Pump Down','Liquid Xenon', '3ms lifetime']
    XRange = [1.0/np.max(multi),np.max(days)]
    YRange = [1E7,1E18]
    Label = ['Time [days]','Total Number of Impurities Left in Plastic']
    dfs.PlotVsTime(XValue=valx[:],YValue=[y1,y2,y4],XRange=XRange,YRange=YRange,Label=Label,Legend=Legend,Log=[False,True])
    plt.savefig(path+'exo_teflon_conc.pdf', dpi=1000, bbox_inches='tight')
    # plt.show()


    C0 = 2.29E5*136.0/32.0       # ppt in mol/mol exo-200 
    C00 = C0 / 136           # ppt in mol/g 
    C00 = C00 * 3.1           # ppt in mol/cm^3
    C00 = C00 * 6.023E23      # ppt in #/cm^3
    C00 = C00 * 1E-12         # units of #/cm^3
    print(C00)

    c1 = y1/6E24*1E12 #impurity concentration ppt [mol/mol]
    c2 = y2/6E24*1E12 #impurity concentration ppt [mol/mol]
    # c3 = y3/6E24*1E12 #impurity concentration ppt [mol/mol]
    c4 = y4/6E24*1E12 #impurity concentration ppt [mol/mol]
    
    YRange = [1E0,1E5]
    Label = ['Time [days]','Impurity Concentration [ppt]']
    dfs.PlotVsTime(XValue=valx,YValue=[c1,c2,c4],XRange=XRange,YRange=YRange,Label=Label,Legend=Legend,Log=[False,True],Loc='lower left')
    plt.savefig(path+'exo_teflon_impur.pdf', dpi=1000, bbox_inches='tight')
    # plt.show()




    c1 = y1/452 # nr of impurities/cm^3 1.4kg/3.1g/mL
    c2 = y2/452 # nr of impurities/cm^3
    # c3 = y3/64500 # nr of impurities/cm^3
    c4 = y4/452 # nr of impurities/cm^3

    Area = 14*14
    f1 = dfs.GetFlowRate(xind[0], Diffs[0], thickness, c1[0], Area)
    for ii,ff in enumerate(f1): 
        if ii>=index1:
            f1[ii] = 0
    f2 = dfs.GetFlowRate(xind[1], Diffs[1], thickness, c2[0], Area)
    for ii,ff in enumerate(f2): 
        if ii>=index2:
            f2[ii] = 0
    # f3 = dfs.GetFlowRate(xind[2], Diffs[2], thickness, c3[0], Area)
    f4 = [1.18E10]*len(xs[0])
    print(C00, c1[-1], c2[-1])
    XRange = [1E-2,np.max(xs)]
    YRange = [1E5,1E20]
    Label = ['Time [days]', r'Outgassing Rate [$\mathrm{Molecules}\,\mathrm{day}^{-1}$]']
    dfs.PlotVsTime(XValue=valx,YValue=[f1,f2,f4],Legend=Legend,Label=Label,XRange=XRange,YRange=YRange,Log=[False,True])
    plt.savefig(path+'exo_teflon_flow.pdf', dpi=1000, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    PlotEXO200Teflon()
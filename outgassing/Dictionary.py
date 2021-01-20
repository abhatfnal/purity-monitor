def GetDictionary():
    #Pump speed is based on https://shop.edwardsvacuum.com/Viewers/Document.ashx?id=2129&lcid=2057
    Dict = {    
        'Helium': {'Mass':4.0, 'Formula': r'$\mathrm{He}$', 'PumpSpeed':78},
        'Carbon': {'Mass':12.0, 'Formula': r'$\mathrm{C}$', 'PumpSpeed':75.5},
        'Methane': {'Mass':16.0, 'Formula': r'$\mathrm{CH}_4$', 'PumpSpeed':75.5},
        'Water': {'Mass':18.0, 'Formula': r'$\mathrm{H}_2\mathrm{O}$', 'PumpSpeed':60},
        'Nitrogen': {'Mass':28.0, 'Formula': r'$\mathrm{N}_2$', 'PumpSpeed':84},
        'Oxygen': {'Mass':32.0, 'Formula': r'$\mathrm{O}_2$', 'PumpSpeed':75.5},
        'Argon': {'Mass':40.0, 'Formula': r'$\mathrm{Ar}$', 'PumpSpeed':80},
        'Ethanol': {'Mass':41.0, 'Formula': r'$\mathrm{C}_2\mathrm{H}_5\mathrm{OH}$', 'PumpSpeed':75.5},
        'Carbondioxide': {'Mass':44.0, 'Formula': r'$\mathrm{CO}_2$', 'PumpSpeed':75.5}}
    return Dict

def GetTurnToLeakRate():
    LeakRate = {
        0.40: {'LeakRate':0.09, 'Error': 0.01},
        1.40: {'LeakRate':0.39, 'Error': 0.04},
        1.60: {'LeakRate':2.36, 'Error': 0.26},
        1.65: {'LeakRate':3.91, 'Error': 0.43},
        1.70: {'LeakRate':7.20, 'Error': 0.79},
        1.75: {'LeakRate':11.73, 'Error': 1.29},
        1.80: {'LeakRate':28.16, 'Error': 3.10},
        1.85: {'LeakRate':64.93, 'Error': 7.14},
        1.90: {'LeakRate':115.78, 'Error': 12.74}}
    return LeakRate
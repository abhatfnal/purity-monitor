import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex' : True,'font.size' : 12,'font.family' : 'lmodern','text.latex.unicode': True,}
plt.rcParams.update(params) 

colors = ['#1f78b4', '#e66101', '#33a02c', '#984ea3']
DriftTime = [76.9776, 104.9937, 161.679636364, 235.649073786, 428.2908] # us
DriftErr = [2.99238765403, 2.82537681558, 4.34909948555, 4.79607649594, 22.6045846146] 
DriftField = [3600.0,1800.0,1040.0,700.0,360.0] # total V across DriftLength

DriftLength = 6.0*25.4 # mm
DriftVelocity = [DriftLength/x for x in DriftTime]
print DriftVelocity

DrfitVelErr = [DriftLength/(x**2) * y for x,y in zip(DriftTime, DriftErr)]
DriftField = [(x*150.0/350.0)/(DriftLength/10.0) for x in DriftField]
print DriftField


fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.ylim(0,500)
plt.xlabel('Drift Field [V/cm]', fontsize=14)
plt.ylabel('Drift Time [$\mu$s]', fontsize=14)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(50))
plt.grid()
plt.errorbar(DriftField, DriftTime, yerr=DriftErr, ms=5, label='', linestyle='none', color=colors[0], marker='o')
plt.legend(loc='upper right')


fig = plt.figure(figsize=(8,6))
ax = fig.gca()
plt.ylim(0,3)
plt.xlabel('Drift Field [V/cm]', fontsize=14)
plt.ylabel('Drift Velocity [mm/$\mu$s]', fontsize=14)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
plt.grid()
plt.errorbar(DriftField, DriftVelocity, yerr=DrfitVelErr, ms=5, label='', linestyle='none', color=colors[0], marker='o')
plt.legend(loc='upper right')

plt.show()
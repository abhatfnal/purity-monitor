import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib, datetime

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex' : True,'font.size' : 12,'font.family' : 'lmodern','text.latex.unicode': True,}
plt.rcParams.update(params) 

dates = ['20190402', '20190608', '20190709', '20190710', '20190802', '20190820', '20190830', '20190923', '20190925', '20191001']
level = [680, 804, 838, 836, 816, 776, 762, 624, 616, 718 ]

dates = [datetime.datetime.strptime(x, '%Y%m%d') for x in dates]

fig = plt.figure(figsize=(12,9))
ax = fig.gca()

ax.grid(b=True, which='major', color='k', linestyle='--', alpha=0.8)
ax.grid(b=True, which='minor', color='grey', linestyle=':')

plt.xlabel('Date' , fontsize=14)
plt.ylabel('Pressure [PSIG]', fontsize=14)

plt.gcf().autofmt_xdate()
formatter = matplotlib.dates.DateFormatter('%b %d %Y')
plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

plt.scatter(dates, level)

datemin = np.datetime64(dates[0], 'm') - np.timedelta64(2, 'm')
datemax = np.datetime64(dates[-1], 'm') + np.timedelta64(1, 'm')

ax.set_xlim(datemin, datemax)

plt.savefig('XenonVsTime.pdf', bbox_inches = 'tight', pad_inches = 0.2)
plt.show()
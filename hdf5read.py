import matplotlib.pyplot as plt
import numpy as np
import h5py


def get_all(name):
    print name

file = '/gpfs/loomis/project/fas/david_moore/aj487/LXeData/20190125/100Vcm_8cm_lxe_anode.h5'
attr = []
f = h5py.File(file, 'r')
attr.append(f.visit(get_all))
print attr
datasetNames = [n for n in f.keys()]
attrNames = [n for n in f.attrs.keys()]
default = f['Default']
attrNames = [n for n in default.attrs.keys()]
print attrNames
print f.items()
for n in datasetNames:
    print(n)
    data = f.get(n)
    print np.array(data)


print
waveformNames = [n for n in f['Default/Waveforms'].keys()]
for n in waveformNames:
    print(n)
    data = np.array(f['Default/Waveforms/%s' % n])
    print len(data)
    print data


    #
    # for key in f.keys():
    #     print key
    #     print np.array(f[key])
        # if(key=='Waveforms'):
        #     kk = np.array(f[key])
        #     for x, i in enumerate(kk):
        #         print np.array(kk[x])

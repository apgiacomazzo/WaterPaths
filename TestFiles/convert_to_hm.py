import numpy as np
from glob import glob

 directories = ['demands', 'inflows']#, 'evaporation']
 for directory in directories:
     files = glob(directory + '/*csv')
     for f in files:
         print 'working on file {}'.format(f)
         data = np.loadtxt(f, delimiter=',') / 1e6
         if data[0][0] > 1000:
             np.savetxt(f, data, delimiter=',', fmt="%.6f")

directories = ['evaporation']
for directory in directories:
    files = glob(directory + '/*csv')
    for f in files:
        print 'working on file {}'.format(f)
        data = np.loadtxt(f, delimiter=',') / 1e2
        if abs(data[0][0]) > 1e-3:
            np.savetxt(f, data, delimiter=',', fmt="%.8f")

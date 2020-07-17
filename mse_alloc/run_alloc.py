from mse_alloc import MseFibAlloc
from matplotlib import pyplot as plt
from astropy.table import Table
import numpy as np

# Without dithering
alloc_nodith = MseFibAlloc(file='cosmo_targets_large_old.csv',
                           doplot=False, meth='fixiter', iternum=1, dither=False, allocfrac=90, spectro='LR')

# Look at results
res = Table.read('TARGETS/results.csv', format='csv')

# Plot Nobsdone as a function of target score (target priority * survey priority)
#plt.scatter(res['priority'] * res['surveypriority'], res['Nobsdone'])
#plt.show()

# Print mean value of Nobsdone for various target score
target_score = np.unique(res['priority'] * res['surveypriority'])
for ts in target_score:
    print(ts, np.mean(res['Nobsdone'][res['priority'] * res['surveypriority'] == ts]))

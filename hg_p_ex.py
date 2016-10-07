import numpy as np
#import scipy as sp
from scipy.stats import hypergeom
import matplotlib as mpl
from matplotlib import pyplot as plt
from collections import OrderedDict

N_stars = 1000.
n_trials = 100.


k_recoverd = np.arange(0,40.5,1)


probs = {}

p = 200.
rv = hypergeom(N_stars, p, n_trials)

for ii, k in enumerate(k_recoverd):
    probs[k] = rv.pmf(k)

P = OrderedDict(sorted(probs.items(), key=lambda t: t[0]))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(P.keys(), P.values())

ax.set_xlabel('True number N of Binaries in 1000 Total stars in Cluster')
ax.set_ylabel('P(N | 25 binaries recoved out of 100 stars')

plt.savefig('Example_BF_hypergeom.png', format='png')


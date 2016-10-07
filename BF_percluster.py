""" BF_percluster

Calculates posterior of P(BF| D,B) for each cluster and makes plot.

"""
import matplotlib as mpl
mpl.use('TkAgg')

from matplotlib import pyplot as plt
import numpy as np
from pandas import read_csv
import os
import bayesian_BF


# Load Data
os.environ['bf_dir'] = '/Users/jquark/projects/karl'
data_dir = os.environ['bf_dir']
sim_results_file = 'insync_simulation_results.txt'
colnames = ['twomassid', 'cluster', 'pdetect', 'binarity', 'baseline']
data = read_csv(os.path.join(data_dir, sim_results_file),
                delim_whitespace=True, header=None, names=colnames)

data.baseline = data.baseline / 86400.  # seconds to days

clusters = data.cluster.unique()
bycluster = data.groupby('cluster')

binary_fraction = {}   # empty dictionaries to hold results by cluster
p_binary_fraction = {}

for ii, cluster in enumerate(clusters):
    obs = bycluster.get_group(cluster)
    BFs, prob_BFs = bayesian_BF.normalized_pBF(obs.pdetect.values,
                                               obs.binarity.values,
                                               BFsteps=1000)
    binary_fraction[cluster] = BFs
    p_binary_fraction[cluster] = prob_BFs


fig = plt.figure()
ax = fig.add_subplot(111)
for cluster in clusters:
    ax.plot(binary_fraction[cluster], p_binary_fraction[cluster],
            label=cluster)
ax.legend()
ax.set_xlabel('Binary Fraction', fontsize=12)
ax.set_ylabel('P(BF| D, B)', fontsize=12)
filename = 'p_BF_percluster'
plt.savefig('{}.png'.format(filename), format='png')

# New figure()
# Pdetect
bins = np.arange(0, 1.01, 0.05)
fig = plt.figure()
ax = fig.add_subplot(111)
for ii, cluster in enumerate(clusters):
    obs = bycluster.get_group(cluster)
    obs['pdetect'].plot(kind='kde', ax=ax, label=cluster)
ax.legend()
filename = 'pdetect_percluster'
plt.savefig('{}.png'.format(filename), format='png')

# New figure()
# Baseline
fig = plt.figure()
ax = fig.add_subplot(111)
for ii, cluster in enumerate(clusters):
    obs = bycluster.get_group(cluster)
    obs['baseline'].plot(kind='hist', cumulative=True, normed=1, ax=ax, label=cluster, histtype='step')
ax.legend()
filename = 'tbaseline_percluster'
plt.savefig('{}.png'.format(filename), format='png')

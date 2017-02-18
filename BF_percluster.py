""" BF_percluster

Calculates posterior of P(BF| D,B) for each cluster and makes plot.

"""
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np
from pandas import read_csv
import bayesian_BF
from os import path


class SimResultsPdetect(object):
    """
    Data Struture Holds Simulation Results
    Reads into pandas DF
    """

    def __init__(self, filename, NDRV=3.0):
        self.filename = filename
        self._ndrv = NDRV  # constant
        colnames = ['twomassid', 'cluster', 'pdetect', 'binarity',
                    'baseline', 'NDRV']
        self.load_data(colnames)
        self.binary_fraction = {}  # Stores results
        self.p_binary_fraction = {}

    def load_data(self, colnames):
        data = read_csv(self.filename,
                        delim_whitespace=True, header=None, names=colnames)
        data.baseline = data.baseline / 86400.  # seconds to days
        self.clusternames = data.cluster.unique()
        bycluster = data.groupby('cluster')
        self._data = data
        self._bycluster = bycluster
        return self

    def cluster_data(self, clustername):
        try:
            return self._bycluster.get_group(clustername)
        except AttributeError:
            self.load_data()
            return self._bycluster.get_group(clustername)

    def compute_BF_prob(self, clustername):
        obs = self.cluster_data(clustername)
        BF_probs = bayesian_BF.normalized_pBF(obs.pdetect.values,
                                              obs.binarity.values,
                                              BFsteps=1000)
        self.binary_fraction[clustername] = BF_probs[0]
        self.p_binary_fraction[clustername] = BF_probs[1]
        return self


class CombinedSimResults(SimResultsPdetect):

    def load_data(self, colnames):
        data = read_csv(self.filename,
                        delim_whitespace=True, header=None, names=colnames)
        binary_flag_file = path.join(self.filename[:self.filename.rfind('/')],
                                     'binarity.{}'.format(int(self._ndrv)))
        binary_flags = np.loadtxt(binary_flag_file)
        # replace bad flags
        data.binarity = binary_flags
        data.baseline = data.baseline / 86400.  # seconds to days
        self.clusternames = data.cluster.unique()
        bycluster = data.groupby('cluster')
        self._data = data
        self._bycluster = bycluster
        return self


def plots_percluster(simresults, plot_dir='plots/'):
    """
    simresults is a class instance of SimResultsPdetect
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for cluster in simresults.clusternames:
        ax.plot(simresults.binary_fraction[cluster],
                simresults.p_binary_fraction[cluster],
                label=cluster)
    ax.legend()
    ax.set_xlabel('Binary Fraction', fontsize=12)
    ax.set_ylabel('P(BF| D, B)', fontsize=12)
    filename = path.join(plot_dir, 'p_BF_percluster')
    plt.savefig('{}_NDRV_{}.png'.format(filename, simresults._ndrv),
                format='png')
    # New figure()
    # Pdetect
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for cluster in simresults.clusternames:
        obs = simresults.cluster_data(cluster)
        obs['pdetect'].plot(kind='kde', ax=ax, label=cluster)
    ax.legend()
    filename = path.join(plot_dir, 'pdetect_percluster')
    plt.savefig('{}_NDRV_{}.png'.format(filename, simresults._ndrv),
                format='png')
    # New figure()
    # Baseline
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for cluster in simresults.clusternames:
        obs = simresults.cluster_data(cluster)
        obs['baseline'].plot(kind='hist', cumulative=True, normed=1, ax=ax,
                             label=cluster, histtype='step')
    ax.legend()
    filename = path.join(plot_dir, 'tbaseline_percluster')
    plt.savefig('{}_NDRV_{}.png'.format(filename, simresults._ndrv),
                format='png')

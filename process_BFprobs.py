import BF_percluster


ndrv3file = './data/insync_simulation_results_ndrv_cutoff_3.txt'
ndrv4file = './data/insync_simulation_results_ndrv_cutoff_4.txt'
ndrv6file = './data/insync_simulation_results_ndrv_cutoff_6.txt'

ndrv3data = BF_percluster.CombinedSimResults(ndrv3file, NDRV=3.0)
ndrv4data = BF_percluster.CombinedSimResults(ndrv4file, NDRV=4.0)
ndrv6data = BF_percluster.CombinedSimResults(ndrv6file, NDRV=6.0)
#ndrv3data = BF_percluster.SimResultsPdetect(ndrv3file, NDRV=3.0)
#ndrv4data = BF_percluster.SimResultsPdetect(ndrv4file, NDRV=4.0)
#ndrv6data = BF_percluster.SimResultsPdetect(ndrv6file, NDRV=6.0)

[ndrv3data.compute_BF_prob(cluster) for cluster in ndrv3data.clusternames]
BF_percluster.plots_percluster(ndrv3data, plot_dir='plots')

[ndrv4data.compute_BF_prob(cluster) for cluster in ndrv4data.clusternames]
BF_percluster.plots_percluster(ndrv4data, plot_dir='plots')

[ndrv6data.compute_BF_prob(cluster) for cluster in ndrv6data.clusternames]
BF_percluster.plots_percluster(ndrv6data, plot_dir='plots')

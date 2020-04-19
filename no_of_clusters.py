# Developed by: Shehu Lab
# Plots the distributions of the number of clusters for the clustering
# techniques

from modules import io
import configparser as cp
import numpy as np

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/kplot.ini")

hierarchical = np.loadtxt(conf['init']['file1path'], usecols=0)
kmeans = np.loadtxt(conf['init']['file2path'], usecols=0)
gmm = np.loadtxt(conf['init']['file3path'], usecols=0)
gromacs = np.loadtxt(conf['init']['file4path'], usecols=0)

# Plot histogram of hierarchical Ks
io.plot_rmsd_distribution(
    hierarchical, 0, 100, 1,
    "Number of Clusters",
    "Frequency", conf['distribution']['color'],
    conf['output']['fileName'] + "hierarchical_K_all"
)

# Plot histogram of kmeans Ks
io.plot_rmsd_distribution(
    kmeans, 0, 100, 1,
    "Number of Clusters",
    "Frequency", conf['distribution']['color'],
    conf['output']['fileName'] + "kmeans_K_all"
)

# Plot histogram of gmm Ks
io.plot_rmsd_distribution(
    gmm, 0, 100, 1,
    "Number of Clusters",
    "Frequency", conf['distribution']['color'],
    conf['output']['fileName'] + "gmm_K_all"
)

# Plot histogram of gromos Ks
io.plot_rmsd_distribution(
    gromacs, 0, 100, 1,
    "Number of Clusters",
    "Frequency", conf['distribution']['color'],
    conf['output']['fileName'] + "gromacs_K_all"
)

# Developed by: Shehu Lab
# Applies GMM clustering and generates reduced ensemble

import configparser as cp
from sklearn.mixture import GaussianMixture
import numpy as np


# ************
# Clustering *
# ************

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/gmm.ini")

# Read data files
data = np.loadtxt(
    conf['init']['dataInputPath'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)

np.seterr(divide='ignore', invalid='ignore')

# Perform clustering
num_of_points = np.size(data, 0)

clusters = GaussianMixture(
        int(conf['init']['numClusters']), covariance_type='full', random_state=42, max_iter=1000
    ).fit_predict(data)

# ************************
# Reduced Pool Selection *
# ************************

# Find the count of each unique label
unique_labels = np.unique(clusters)

# Read energy and rmsd of each point
energies = np.loadtxt(conf['init']['dataInputPath'], usecols=0)
rmsds = np.loadtxt(conf['init']['dataInputPath'], usecols=1)

# Extract points from each cluster based on energy score
selected = []

for i in unique_labels:
    indices = np.where(clusters == i)[0]
    group_energies = energies[indices]

    unique_energies, unique_indices = np.unique(
        group_energies, return_index=True
    )

    unique_energies = np.around(unique_energies, decimals=2)

    selection_indices = indices[unique_indices]

    taken = False
    energy_level = np.min(unique_energies)

    for j in range(len(unique_energies)):
        energy = unique_energies[j]

        if energy > energy_level:
            energy_level = energy
            taken = False

        if energy == energy_level:
            if not taken:
                selected.append(selection_indices[j])
                taken = True

# Write output
output_file = conf['output']['outputFile']

f = open(output_file, "a")

for i in selected:
    f.write("%s %s\n" % (energies[i], rmsds[i]))

f.close()

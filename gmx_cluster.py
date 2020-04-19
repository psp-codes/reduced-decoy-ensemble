# Developed by: Shehu Lab
# Applies gmx-cluster clustering and generates reduced ensemble

import configparser as cp
import numpy as np
from scipy.spatial import distance
import math

# ****************
# Initialization *
# ****************

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/gmx_cluster.ini")

# Read data files
usr_features = np.loadtxt(
    conf['init']['dataInputPath'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)

np.seterr(divide='ignore', invalid='ignore')

# ************
# Clustering *
# ************

# Find pairwise distance
distance_pairs = []
max_distance = -math.inf
min_distance = math.inf

no_of_points = np.size(usr_features, 0)

for i in range(no_of_points):
    distance_individual = []
    for j in range(no_of_points):
        dist = distance.euclidean(usr_features[i], usr_features[j])

        if dist > max_distance:
            max_distance = dist

        distance_individual.append(dist)

    distance_pairs.append(distance_individual)

dist_pairs = np.array(distance_pairs)

dist_pairs = dist_pairs / max_distance

cutoff = float(conf['init']['cutoff'])

clusters = []
clustered = []
indices = [*range(no_of_points)]

while True:

    max_neighbor_count = 0
    max_neighbors = []

    for i in indices:
        neighbors = []
        neighbor_count = 0

        for j in indices:  # change to indices
            if dist_pairs[i][j] <= cutoff:
                neighbor_count += 1
                neighbors.append(j)

        if neighbor_count > max_neighbor_count:
            max_neighbor_count = neighbor_count
            max_neighbors = list(neighbors)

    clusters.append(max_neighbors)

    for i in max_neighbors:
        indices.remove(i)

    if len(indices) == 0:
        break

# ************************
# Reduced Pool Selection *
# ************************

# Read energy and rmsd of each point
energies = np.loadtxt(conf['init']['dataInputPath'], usecols=0)
rmsds = np.loadtxt(conf['init']['dataInputPath'], usecols=1)

# Extract points from each cluster based on energy score
selected = []

for cluster in clusters:
    cluster_array = np.array(cluster)

    group_energies = energies[cluster_array]

    unique_energies, unique_indices = np.unique(
        group_energies, return_index=True
    )

    unique_energies = np.around(unique_energies, decimals=2)

    selection_indices = cluster_array[unique_indices]

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




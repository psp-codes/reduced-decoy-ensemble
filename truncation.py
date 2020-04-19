# Developed by: Shehu Lab
# Applies truncation and generates reduced ensemble

import configparser as cp
import numpy as np


# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/truncate.ini")

# ************************
# Reduced Pool Selection *
# ************************

energies = np.loadtxt(conf['init']['dataInputPath'], usecols=0)
rmsds = np.loadtxt(conf['init']['dataInputPath'], usecols=1)

reduced_hier = np.loadtxt(conf['init']['file1path'], usecols=0)
reduced_kmeans = np.loadtxt(conf['init']['file2path'], usecols=0)
reduced_gmm = np.loadtxt(conf['init']['file3path'], usecols=0)
reduced_gromacs = np.loadtxt(conf['init']['file4path'], usecols=0)

lengths = [reduced_hier.size, reduced_kmeans.size, reduced_gmm.size,
           reduced_gromacs.size]
output_length = min(lengths)

# Extract points based on energy score
sorted_energies = np.argsort(energies)

# Write output
output_file = conf['output']['outputFile']

f = open(output_file, "a")

for i in range(output_length):
    f.write("%s %s\n" % (energies[sorted_energies[i]], rmsds[sorted_energies[i]]))

f.close()

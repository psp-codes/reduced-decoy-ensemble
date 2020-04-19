# Developed by: Shehu Lab
# Calculates and plots energy comparison results

from modules import io
import configparser as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyrosetta as pr

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/energy_comparison.ini")

# Read data files
full_energy1 = np.loadtxt(conf['full_size']['file1path'], usecols=0)
full_energy2 = np.loadtxt(conf['full_size']['file2path'], usecols=0)
full_energy3 = np.loadtxt(conf['full_size']['file3path'], usecols=0)
full_energy4 = np.loadtxt(conf['full_size']['file4path'], usecols=0)
full_energy5 = np.loadtxt(conf['full_size']['file5path'], usecols=0)

reduced_energy1 = np.loadtxt(conf['reduced']['file1path'], usecols=0)
reduced_energy2 = np.loadtxt(conf['reduced']['file2path'], usecols=0)
reduced_energy3 = np.loadtxt(conf['reduced']['file3path'], usecols=0)
reduced_energy4 = np.loadtxt(conf['reduced']['file4path'], usecols=0)
reduced_energy5 = np.loadtxt(conf['reduced']['file5path'], usecols=0)

reduced_energy1_hier = np.loadtxt(conf['reduced_hier']['file1path'], usecols=0)
reduced_energy2_hier = np.loadtxt(conf['reduced_hier']['file2path'], usecols=0)
reduced_energy3_hier = np.loadtxt(conf['reduced_hier']['file3path'], usecols=0)
reduced_energy4_hier = np.loadtxt(conf['reduced_hier']['file4path'], usecols=0)
reduced_energy5_hier = np.loadtxt(conf['reduced_hier']['file5path'], usecols=0)

reduced_energy1_gmm = np.loadtxt(conf['reduced_gmm']['file1path'], usecols=0)
reduced_energy2_gmm = np.loadtxt(conf['reduced_gmm']['file2path'], usecols=0)
reduced_energy3_gmm = np.loadtxt(conf['reduced_gmm']['file3path'], usecols=0)
reduced_energy4_gmm = np.loadtxt(conf['reduced_gmm']['file4path'], usecols=0)
reduced_energy5_gmm = np.loadtxt(conf['reduced_gmm']['file5path'], usecols=0)

reduced_energy1_gromacs = np.loadtxt(conf['reduced_gromacs']['file1path'], usecols=0)
reduced_energy2_gromacs = np.loadtxt(conf['reduced_gromacs']['file2path'], usecols=0)
reduced_energy3_gromacs = np.loadtxt(conf['reduced_gromacs']['file3path'], usecols=0)
reduced_energy4_gromacs = np.loadtxt(conf['reduced_gromacs']['file4path'], usecols=0)
reduced_energy5_gromacs = np.loadtxt(conf['reduced_gromacs']['file5path'], usecols=0)

reduced_energy1_truncation = np.loadtxt(conf['reduced_truncation']['file1path'], usecols=0)
reduced_energy2_truncation = np.loadtxt(conf['reduced_truncation']['file2path'], usecols=0)
reduced_energy3_truncation = np.loadtxt(conf['reduced_truncation']['file3path'], usecols=0)
reduced_energy4_truncation = np.loadtxt(conf['reduced_truncation']['file4path'], usecols=0)
reduced_energy5_truncation = np.loadtxt(conf['reduced_truncation']['file5path'], usecols=0)

full_energies = np.concatenate(
    (full_energy1, full_energy2, full_energy3, full_energy4, full_energy5)
)

reduced_energies_km = np.concatenate(
    (reduced_energy1, reduced_energy2, reduced_energy3, reduced_energy4,
     reduced_energy5)
)

reduced_energies_hier = np.concatenate(
    (reduced_energy1_hier, reduced_energy2_hier, reduced_energy3_hier, reduced_energy4_hier, reduced_energy5_hier)
)

reduced_energies_gmm = np.concatenate(
    (reduced_energy1_gmm, reduced_energy2_gmm, reduced_energy3_gmm, reduced_energy4_gmm, reduced_energy5_gmm)
)

reduced_energies_gromacs = np.concatenate(
    (reduced_energy1_gromacs, reduced_energy2_gromacs, reduced_energy3_gromacs, reduced_energy4_gromacs, reduced_energy5_gromacs)
)

reduced_energies_truncation = np.concatenate(
    (reduced_energy1_truncation, reduced_energy2_truncation, reduced_energy3_truncation, reduced_energy4_truncation, reduced_energy5_truncation)
)

# Outputs

print("Min energy original: %.2f" % full_energies.min())
print("Min energy kmeans: %.2f" % reduced_energies_km.min())
print("Min energy gmm: %.2f" % reduced_energies_gmm.min())
print("Min energy hier: %.2f" % reduced_energies_hier.min())
print("Min energy gromacs: %.2f" % reduced_energies_gromacs.min())
print("Min energy truncation: %.2f" % reduced_energies_truncation.min())

print("")

print("Avg energy original: %.2f" % full_energies.mean())
print("Avg energy kmeans: %.2f" % reduced_energies_km.mean())
print("Avg energy gmm: %.2f" % reduced_energies_gmm.mean())
print("Avg energy hier: %.2f" % reduced_energies_hier.mean())
print("Avg energy gromacs: %.2f" % reduced_energies_gromacs.mean())
print("Avg energy truncation: %.2f" % reduced_energies_truncation.mean())

print("")

print("Diameter energy original: %.2f" % (full_energies.max() - full_energies.min()))
print("Diameter energy kmeans: %.2f" % (reduced_energies_km.max() - reduced_energies_km.min()))
print("Diameter energy gmm: %.2f" % (reduced_energies_gmm.max() - reduced_energies_gmm.min()))
print("Diameter energy hier: %.2f" % (reduced_energies_hier.max() - reduced_energies_hier.min()))
print("Diameter energy gromacs: %.2f" % (reduced_energies_gromacs.max() - reduced_energies_gromacs.min()))
print("Diameter energy truncation: %.2f" % (reduced_energies_truncation.max() - reduced_energies_truncation.min()))


# Plot energy distribution of all pools

mpl.rcParams['mathtext.default'] = 'regular'

output_filename = conf['output']['energyDist'] + "_energy_dist"

bin_array = np.arange(full_energies.min()-1, full_energies.max() + 5, 2)

y1, bin_edges = np.histogram(full_energies, bins=bin_array)
y2, bin_edges = np.histogram(reduced_energies_km, bins=bin_array)
y3, bin_edges = np.histogram(reduced_energies_hier, bins=bin_array)
y4, bin_edges = np.histogram(reduced_energies_gmm, bins=bin_array)
y5, bin_edges = np.histogram(reduced_energies_gromacs, bins=bin_array)

y1 = y1 / (np.max(y1) - np.min(y1))
y2 = y2 / (np.max(y2) - np.min(y2))
y3 = y3 / (np.max(y3) - np.min(y3))
y4 = y4 / (np.max(y4) - np.min(y4))
y5 = y5 / (np.max(y5) - np.min(y5))

plt.figure(figsize=(8, 5))

plt.plot(bin_edges[:-1], y1, linewidth=2, color="#d53e4f", linestyle="-")
plt.plot(bin_edges[:-1], y3, linewidth=2, color="#1a9850", linestyle=(0, (1, 2)))
plt.plot(bin_edges[:-1], y2, linewidth=2, color="#542788", linestyle=(0, (5, 4)))
plt.plot(bin_edges[:-1], y4, linewidth=2, color="#a6611a", linestyle=(0, (1, 1)))
plt.plot(bin_edges[:-1], y5, linewidth=2, color="#3288bd", linestyle=(0, (5, 2)))

plt.legend(["" + r'$\Omega_{gen}$', r'$\Omega_{red}$' + " (hierarchical)",
            r'$\Omega_{red}$' + " (k-means)", '$\Omega_{red}$' + " (GMM)",
            '$\Omega_{red}$' + " (gmx-cluster-usr)"], bbox_to_anchor=(0.74, 0.308))
plt.xlabel('Rosetta score4 Energy (REU)')
plt.ylabel("Frequency")

# plt.show()

plt.savefig(output_filename + ".png", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)
plt.savefig(output_filename + ".pdf", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)

# Developed by: Shehu Lab
# Calculates and plots lRMSD distance results
# Plots visualization of the ensembles

from modules import io
import configparser as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/evaluate.ini")

# Read data files
full_energy1 = np.loadtxt(conf['full_size']['file1path'], usecols=0)
full_rmsd1 = np.loadtxt(conf['full_size']['file1path'], usecols=1)
full_energy2 = np.loadtxt(conf['full_size']['file2path'], usecols=0)
full_rmsd2 = np.loadtxt(conf['full_size']['file2path'], usecols=1)
full_energy3 = np.loadtxt(conf['full_size']['file3path'], usecols=0)
full_rmsd3 = np.loadtxt(conf['full_size']['file3path'], usecols=1)
full_energy4 = np.loadtxt(conf['full_size']['file4path'], usecols=0)
full_rmsd4 = np.loadtxt(conf['full_size']['file4path'], usecols=1)
full_energy5 = np.loadtxt(conf['full_size']['file5path'], usecols=0)
full_rmsd5 = np.loadtxt(conf['full_size']['file5path'], usecols=1)

reduced_energy1 = np.loadtxt(conf['reduced_kmeans']['file1path'], usecols=0)
reduced_rmsd1 = np.loadtxt(conf['reduced_kmeans']['file1path'], usecols=1)
reduced_energy2 = np.loadtxt(conf['reduced_kmeans']['file2path'], usecols=0)
reduced_rmsd2 = np.loadtxt(conf['reduced_kmeans']['file2path'], usecols=1)
reduced_energy3 = np.loadtxt(conf['reduced_kmeans']['file3path'], usecols=0)
reduced_rmsd3 = np.loadtxt(conf['reduced_kmeans']['file3path'], usecols=1)
reduced_energy4 = np.loadtxt(conf['reduced_kmeans']['file4path'], usecols=0)
reduced_rmsd4 = np.loadtxt(conf['reduced_kmeans']['file4path'], usecols=1)
reduced_energy5 = np.loadtxt(conf['reduced_kmeans']['file5path'], usecols=0)
reduced_rmsd5 = np.loadtxt(conf['reduced_kmeans']['file5path'], usecols=1)

reduced_rmsd1_hier = np.loadtxt(conf['reduced_hier']['file1path'], usecols=1)
reduced_rmsd2_hier = np.loadtxt(conf['reduced_hier']['file2path'], usecols=1)
reduced_rmsd3_hier = np.loadtxt(conf['reduced_hier']['file3path'], usecols=1)
reduced_rmsd4_hier = np.loadtxt(conf['reduced_hier']['file4path'], usecols=1)
reduced_rmsd5_hier = np.loadtxt(conf['reduced_hier']['file5path'], usecols=1)

reduced_energy1_hier = np.loadtxt(conf['reduced_hier']['file1path'], usecols=0)
reduced_energy2_hier = np.loadtxt(conf['reduced_hier']['file2path'], usecols=0)
reduced_energy3_hier = np.loadtxt(conf['reduced_hier']['file3path'], usecols=0)
reduced_energy4_hier = np.loadtxt(conf['reduced_hier']['file4path'], usecols=0)
reduced_energy5_hier = np.loadtxt(conf['reduced_hier']['file5path'], usecols=0)

reduced_rmsd1_gmm = np.loadtxt(conf['reduced_gmm']['file1path'], usecols=1)
reduced_rmsd2_gmm = np.loadtxt(conf['reduced_gmm']['file2path'], usecols=1)
reduced_rmsd3_gmm = np.loadtxt(conf['reduced_gmm']['file3path'], usecols=1)
reduced_rmsd4_gmm = np.loadtxt(conf['reduced_gmm']['file4path'], usecols=1)
reduced_rmsd5_gmm = np.loadtxt(conf['reduced_gmm']['file5path'], usecols=1)

reduced_energy1_gmm = np.loadtxt(conf['reduced_gmm']['file1path'], usecols=0)
reduced_energy2_gmm = np.loadtxt(conf['reduced_gmm']['file2path'], usecols=0)
reduced_energy3_gmm = np.loadtxt(conf['reduced_gmm']['file3path'], usecols=0)
reduced_energy4_gmm = np.loadtxt(conf['reduced_gmm']['file4path'], usecols=0)
reduced_energy5_gmm = np.loadtxt(conf['reduced_gmm']['file5path'], usecols=0)

reduced_rmsd1_gromacs = np.loadtxt(conf['reduced_gmx']['file1path'], usecols=1)
reduced_rmsd2_gromacs = np.loadtxt(conf['reduced_gmx']['file2path'], usecols=1)
reduced_rmsd3_gromacs = np.loadtxt(conf['reduced_gmx']['file3path'], usecols=1)
reduced_rmsd4_gromacs = np.loadtxt(conf['reduced_gmx']['file4path'], usecols=1)
reduced_rmsd5_gromacs = np.loadtxt(conf['reduced_gmx']['file5path'], usecols=1)

reduced_energy1_gromacs = np.loadtxt(conf['reduced_gmx']['file1path'], usecols=0)
reduced_energy2_gromacs = np.loadtxt(conf['reduced_gmx']['file2path'], usecols=0)
reduced_energy3_gromacs = np.loadtxt(conf['reduced_gmx']['file3path'], usecols=0)
reduced_energy4_gromacs = np.loadtxt(conf['reduced_gmx']['file4path'], usecols=0)
reduced_energy5_gromacs = np.loadtxt(conf['reduced_gmx']['file5path'], usecols=0)

reduced_energy1_truncation = np.loadtxt(conf['reduced_truncation']['file1path'], usecols=0)
reduced_energy2_truncation = np.loadtxt(conf['reduced_truncation']['file2path'], usecols=0)
reduced_energy3_truncation = np.loadtxt(conf['reduced_truncation']['file3path'], usecols=0)
reduced_energy4_truncation = np.loadtxt(conf['reduced_truncation']['file4path'], usecols=0)
reduced_energy5_truncation = np.loadtxt(conf['reduced_truncation']['file5path'], usecols=0)

reduced_rmsd1_truncation = np.loadtxt(conf['reduced_truncation']['file1path'], usecols=1)
reduced_rmsd2_truncation = np.loadtxt(conf['reduced_truncation']['file2path'], usecols=1)
reduced_rmsd3_truncation = np.loadtxt(conf['reduced_truncation']['file3path'], usecols=1)
reduced_rmsd4_truncation = np.loadtxt(conf['reduced_truncation']['file4path'], usecols=1)
reduced_rmsd5_truncation = np.loadtxt(conf['reduced_truncation']['file5path'], usecols=1)

full_rmsds = np.concatenate(
    (full_rmsd1, full_rmsd2, full_rmsd3, full_rmsd4, full_rmsd5)
)
full_energies = np.concatenate(
    (full_energy1, full_energy2, full_energy3, full_energy4, full_energy5)
)

reduced_rmsds = np.concatenate(
    (reduced_rmsd1, reduced_rmsd2, reduced_rmsd3, reduced_rmsd4, reduced_rmsd5)
)
reduced_energies = np.concatenate(
    (reduced_energy1, reduced_energy2, reduced_energy3, reduced_energy4,
     reduced_energy5)
)

reduced_rmsds_hier = np.concatenate(
    (reduced_rmsd1_hier, reduced_rmsd2_hier, reduced_rmsd3_hier, reduced_rmsd4_hier, reduced_rmsd5_hier)
)

reduced_energies_hier = np.concatenate(
    (reduced_energy1_hier, reduced_energy2_hier, reduced_energy3_hier, reduced_energy4_hier, reduced_energy5_hier)
)

reduced_rmsds_gmm = np.concatenate(
    (reduced_rmsd1_gmm, reduced_rmsd2_gmm, reduced_rmsd3_gmm, reduced_rmsd4_gmm, reduced_rmsd5_gmm)
)

reduced_energies_gmm = np.concatenate(
    (reduced_energy1_gmm, reduced_energy2_gmm, reduced_energy3_gmm, reduced_energy4_gmm, reduced_energy5_gmm)
)

reduced_rmsds_gromacs = np.concatenate(
    (reduced_rmsd1_gromacs, reduced_rmsd2_gromacs, reduced_rmsd3_gromacs, reduced_rmsd4_gromacs, reduced_rmsd5_gromacs)
)

reduced_energies_gromacs = np.concatenate(
    (reduced_energy1_gromacs, reduced_energy2_gromacs, reduced_energy3_gromacs, reduced_energy4_gromacs, reduced_energy5_gromacs)
)

reduced_energies_truncation = np.concatenate(
    (reduced_energy1_truncation, reduced_energy2_truncation, reduced_energy3_truncation, reduced_energy4_truncation, reduced_energy5_truncation)
)

reduced_rmsds_truncation = np.concatenate(
    (reduced_rmsd1_truncation, reduced_rmsd2_truncation, reduced_rmsd3_truncation, reduced_rmsd4_truncation, reduced_rmsd5_truncation)
)

# Outputs

print("Size Hier", reduced_rmsds_hier.size)
print("Size Kmeans", reduced_rmsds.size)
print("Size GMM", reduced_rmsds_gmm.size)
print("Size Gromos", reduced_rmsds_gromacs.size)

print("Min lRMSD Full %.2f" % full_rmsds.min())
print("Min lRMSD Hier %.2f" % reduced_rmsds_hier.min())
print("Min lRMSD Kmeans %.2f" % reduced_rmsds.min())
print("Min lRMSD GMM %.2f" % reduced_rmsds_gmm.min())
print("Min lRMSD Gromos %.2f" % reduced_rmsds_gromacs.min())
print("Min lRMSD Truncation %.2f" % reduced_rmsds_truncation.min())

print("Mean lRMSD Full %.2f" % full_rmsds.mean())
print("Mean lRMSD Hier %.2f" % reduced_rmsds_hier.mean())
print("Mean lRMSD Kmeans %.2f" % reduced_rmsds.mean())
print("Mean lRMSD GMM %.2f" % reduced_rmsds_gmm.mean())
print("Mean lRMSD Gromos %.2f" % reduced_rmsds_gromacs.mean())

print("Stddev lRMSD Full %.2f" % full_rmsds.std())
print("Stddev lRMSD Hier %.2f" % reduced_rmsds_hier.std())
print("Stddev lRMSD Kmeans %.2f" % reduced_rmsds.std())
print("Stddev lRMSD GMM %.2f" % reduced_rmsds_gmm.std())
print("Stddev lRMSD Gromos %.2f" % reduced_rmsds_gromacs.std())

# Plot landscape
max_value = int(np.max(full_rmsds))
mpl.rcParams['mathtext.default'] = 'regular'

io.plot_energy_vs_rmsd(
    full_rmsds, full_energies, reduced_rmsds, reduced_energies,
    range(0, max_value, 2),
    'C' + r'$\alpha$' + ' lRMSD to native structure (' + r'$\AA$' + ')',
    "Rosetta score4 energy (REU)",
    float(conf['energyVsRMSD']['pointSize']),
    conf['energyVsRMSD']['set1color'],
    conf['energyVsRMSD']['set2color'],
    conf['energyVsRMSD']['set1marker'],
    conf['energyVsRMSD']['set2marker'],
    conf['output']['fileName'] + "_landscape_km"
)

io.plot_energy_vs_rmsd(
    full_rmsds, full_energies, reduced_rmsds_hier, reduced_energies_hier,
    range(0, max_value, 2),
    'C' + r'$\alpha$' + ' lRMSD to native structure (' + r'$\AA$' + ')',
    "Rosetta score4 energy (REU)",
    float(conf['energyVsRMSD']['pointSize']),
    conf['energyVsRMSD']['set1color'],
    conf['energyVsRMSD']['set2color'],
    conf['energyVsRMSD']['set1marker'],
    conf['energyVsRMSD']['set2marker'],
    conf['output']['fileName'] + "_landscape_hier"
)

io.plot_energy_vs_rmsd(
    full_rmsds, full_energies, reduced_rmsds_gmm, reduced_energies_gmm,
    range(0, max_value, 2),
    'C' + r'$\alpha$' + ' lRMSD to native structure (' + r'$\AA$' + ')',
    "Rosetta score4 energy (REU)",
    float(conf['energyVsRMSD']['pointSize']),
    conf['energyVsRMSD']['set1color'],
    conf['energyVsRMSD']['set2color'],
    conf['energyVsRMSD']['set1marker'],
    conf['energyVsRMSD']['set2marker'],
    conf['output']['fileName'] + "_landscape_gmm"
)

io.plot_energy_vs_rmsd(
    full_rmsds, full_energies, reduced_rmsds_gromacs, reduced_energies_gromacs,
    range(0, max_value, 2),
    'C' + r'$\alpha$' + ' lRMSD to native structure (' + r'$\AA$' + ')',
    "Rosetta score4 energy (REU)",
    float(conf['energyVsRMSD']['pointSize']),
    conf['energyVsRMSD']['set1color'],
    conf['energyVsRMSD']['set2color'],
    conf['energyVsRMSD']['set1marker'],
    conf['energyVsRMSD']['set2marker'],
    conf['output']['fileName'] + "_landscape_gromos"
)

# Plot lRMSD distribution of all pools

output_filename = conf['output']['rmsdDist'] + "_rmsd_dist"

bin_array = np.arange(0, 25 + 1, 1)

y1, bin_edges = np.histogram(full_rmsds, bins=bin_array)
y2, bin_edges = np.histogram(reduced_rmsds, bins=bin_array)
y3, bin_edges = np.histogram(reduced_rmsds_hier, bins=bin_array)
y4, bin_edges = np.histogram(reduced_rmsds_gmm, bins=bin_array)
y5, bin_edges = np.histogram(reduced_rmsds_gromacs, bins=bin_array)

y1 = y1 / (np.max(y1) - np.min(y1))
y2 = y2 / (np.max(y2) - np.min(y2))
y3 = y3 / (np.max(y3) - np.min(y3))
y4 = y4 / (np.max(y4) - np.min(y4))
y5 = y5 / (np.max(y5) - np.min(y5))

plt.plot(bin_edges[:-1], y1, linewidth=2, color="#d53e4f", linestyle="-")
plt.plot(bin_edges[:-1], y3, linewidth=2, color="#1a9850", linestyle=(0, (1, 2)))
plt.plot(bin_edges[:-1], y2, linewidth=2, color="#542788", linestyle=(0, (5, 4)))
plt.plot(bin_edges[:-1], y4, linewidth=2, color="#a6611a", linestyle=(0, (1, 1)))
plt.plot(bin_edges[:-1], y5, linewidth=2, color="#3288bd", linestyle=(0, (5, 2)))

plt.legend(["Original pool", "Reduced pool Hierarchical",
            "Reduced pool K-means", "Reduced pool GMM",
            "Reduced pool Gmx-cluster"])
plt.xlabel('C' + r'$\alpha$' + ' lRMSD to native structure (' + r'$\AA$' + ')')
plt.ylabel("Frequency")

# plt.show()

plt.savefig(output_filename + ".png", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)
plt.savefig(output_filename + ".pdf", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)

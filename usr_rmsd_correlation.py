# Developed by: Shehu Lab
# Calculates USR scores and shows correlation between USR scores
# and lRMSDs to the native structure of the decoys

from modules import feature_extraction as fe
import pyrosetta as pr
import configparser as cp
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import scipy
import matplotlib as mpl

pr.init()

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/correlation.ini")

# Read data file
rmsds1 = np.loadtxt(conf['init']['dataInputPath1'], usecols=1)
rmsds2 = np.loadtxt(conf['init']['dataInputPath2'], usecols=1)
rmsds3 = np.loadtxt(conf['init']['dataInputPath3'], usecols=1)
rmsds4 = np.loadtxt(conf['init']['dataInputPath4'], usecols=1)
rmsds5 = np.loadtxt(conf['init']['dataInputPath5'], usecols=1)

# Native pose
native_pose = pr.pose_from_pdb(conf['init']['pdbPath'])

native_usr = np.array(fe.similarity_metrics(native_pose))

usr_features1 = np.loadtxt(
    conf['init']['dataInputPath1'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)
usr_features2 = np.loadtxt(
    conf['init']['dataInputPath2'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)
usr_features3 = np.loadtxt(
    conf['init']['dataInputPath3'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)
usr_features4 = np.loadtxt(
    conf['init']['dataInputPath4'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)
usr_features5 = np.loadtxt(
    conf['init']['dataInputPath5'],
    usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
)

rmsds = np.concatenate(
    (rmsds1, rmsds2, rmsds3, rmsds4, rmsds5)
)

usr_features = np.concatenate(
    (usr_features1, usr_features2, usr_features3, usr_features4, usr_features5)
)

# Calculate USR score
usrs = []
for usr in usr_features:
    usrs.append(distance.euclidean(usr, native_usr))

usr_scores = np.array(usrs)

# Find correlation
pearson = np.corrcoef(usr_scores, rmsds)
print("Pearson: {:.2f}".format(pearson[0][1]))

# Plot lRMSD vs. USR
mpl.rcParams['mathtext.default'] = 'regular'

slope, intercept, r, p, stderr = scipy.stats.linregress(usr_scores, rmsds)

plt.plot(usr_scores, rmsds, linewidth=0, marker='o', markersize=2)
plt.plot(usr_scores, intercept + slope * usr_scores, linewidth=1)
plt.xlabel('USR Score', fontsize=13)
plt.ylabel('CA' + ' lRMSD to Native Structure (' + r'$\AA$' + ')', fontsize=13)

plt.xticks(range(0, 900, 100))
plt.yticks(range(0, 35, 5))
plt.text(0, 28, "Pearson’s coefficient = {:.2f}".format(pearson[0][1]),
         fontsize=16)
plt.savefig(conf['output']['outputFile'] + ".pdf", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)

plt.savefig(conf['output']['outputFile'] + ".png", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)

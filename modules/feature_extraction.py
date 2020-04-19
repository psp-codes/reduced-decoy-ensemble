# feature_extraction.py
# since: 11/2018
# Developed by: Shehu Lab

"""Module for extracting features of a conformation.

This module provides methods to extract the USR features of a structure.

Available Functions:
- similarity_metrics: Returns USR moments of ctd, fct, and ftf metrics
    of a conformation.
"""

import math
import pyrosetta.toolbox.extract_coords_pose as pte
import numpy as np
import scipy.stats as st


def similarity_metrics(pose):
    """Returns moments of ctd, cst, fct, and ftf metrics of a
    conformation.

    This function calculates Ultrafast Shape Recognition Matrics (i.e.
    moments of ctd, cst, fct, and ftf distances) of a given structure.

    Args:
        pose: A pyrosetta Pose object containing a structure.

    Returns:
        A list containing 12 moments of the USR distances.
        The order is:
        [ctd_distances_mean, ctd_distances_variance,
        ctd_distances skewness, cst_distances_mean, ....]
    """
    xyz = pte.extract_coordinates_from_pose_1x3(pose)
    no_of_atoms = len(xyz)

    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0

    for i in xyz:
        sum_x += i[0]
        sum_y += i[1]
        sum_z += i[2]

    # Find molecular centroid
    center_x = sum_x / no_of_atoms
    center_y = sum_y / no_of_atoms
    center_z = sum_z / no_of_atoms

    # Calculate moments of ctd distances
    max_distance = 0
    min_distance = math.inf
    fct = []
    cst = []
    distances_list = []

    for i in xyz:
        distance = math.sqrt(
            ((center_x - i[0]) ** 2) + ((center_y - i[1]) ** 2) +
            ((center_z - i[2]) ** 2))
        if distance > max_distance:
            max_distance = distance
            fct = i
        if distance < min_distance:
            min_distance = distance
            cst = i
        distances_list.append(distance)

    distances = np.array(distances_list)
    mean_ctd_distance = np.mean(distances)
    variance_ctd_distance = np.var(distances)
    skewness_ctd_distance = st.skew(distances)

    # Calculate moments of cst distances
    distances_list = []

    for i in xyz:
        distances_list.append(math.sqrt(((cst[0] - i[0]) ** 2)
                                        + ((cst[1] - i[1]) ** 2)
                                        + ((cst[2] - i[2]) ** 2)))

    distances = np.array(distances_list)
    mean_cst_distance = np.mean(distances)
    variance_cst_distance = np.var(distances)
    skewness_cst_distance = st.skew(distances)

    # Calculate moments of fct distances
    max_distance = 0
    ftf = []
    distances_list = []

    for i in xyz:
        distance = math.sqrt(((fct[0] - i[0]) ** 2) + ((fct[1] - i[1]) ** 2) +
                             ((fct[2] - i[2]) ** 2))
        if distance > max_distance:
            max_distance = distance
            ftf = i

        distances_list.append(distance)

    distances = np.array(distances_list)
    mean_fct_distance = np.mean(distances)
    variance_fct_distance = np.var(distances)
    skewness_fct_distance = st.skew(distances)

    # Calculate moments of ftf distances
    distances_list = []

    for i in xyz:
        distances_list.append(math.sqrt(((ftf[0] - i[0]) ** 2) +
                                        ((ftf[1] - i[1]) ** 2) +
                                        ((ftf[2] - i[2]) ** 2)))

    distances = np.array(distances_list)
    mean_ftf_distance = np.mean(distances)
    variance_ftf_distance = np.var(distances)
    skewness_ftf_distance = st.skew(distances)

    return [mean_ctd_distance, variance_ctd_distance, skewness_ctd_distance,
            mean_cst_distance, variance_cst_distance, skewness_cst_distance,
            mean_fct_distance, variance_fct_distance, skewness_fct_distance,
            mean_ftf_distance, variance_ftf_distance, skewness_ftf_distance]

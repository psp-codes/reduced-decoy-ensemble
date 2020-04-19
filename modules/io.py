# io.py
# since: 10/2018
# Developed by: Shehu Lab

"""Module for the input/output operations.

Contains input/output utility functions needed for the project.

Available functions:
- get_sequence_from_fasta: Returns the sequence of a protein from it's
    fasta file.
- plot_rmsd_distribution: Plots frequencies for Ca-lRMSDs to the native
    structure.
- plot_energy_vs_rmsd: Plots energies for Ca-lRMSDs to the native
    structure on two sets of decoys.
"""

import matplotlib.pyplot as plt
import numpy as np


def get_sequence_from_fasta(fasta_file):
    """Returns the sequence of a protein from its fasta file

    Args:
        fasta_file: A string containing the fasta file path with the
        file name and extension.

    Raises:
        ValueError: if fasta_file path is empty.

    Returns:
        A string containing the sequence of the protein.
    """

    if not fasta_file:
        raise ValueError("Path cannot be empty.")

    handle = open(fasta_file, 'r')
    try:
        sequence = handle.readlines()
    finally:
        handle.close()

    sequence = [line.strip() for line in sequence if not '>' in line]
    sequence = ''.join(sequence)

    return sequence


def plot_rmsd_distribution(rmsds, x_min, x_max, x_increment, x_label, y_label,
                           plot_color, output_filename):
    """Plots frequencies for Ca-lRMSDs to the native structure. Saves the
    plot on .png and .pdf file formats.

    Args:
        rmsds: A numpy array containing all the Ca-RMSD values.
        x_min: A number indicating the minimum value for X axis.
        x_max: An number indicating the maximum value for X axis.
        x_increment: A float indicating the bin length.
        x_label: A string indicating the label along the X axis.
        y_label: A string indicating the label along the Y axis.
        plot_color: A string indicating the color of the plot.
        output_filename: A string containing the name of the output
            file.

    Returns:
        None
    """

    # Create bin array
    bin_array = np.arange(x_min, x_max + x_increment, x_increment)

    # Create the histogram
    y, bin_edges = np.histogram(rmsds, bins=bin_array)

    # Plot
    plt.plot(bin_edges[:-1], y, color=plot_color)

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.savefig(output_filename + ".pdf", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)
    plt.savefig(output_filename + ".png", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)


def plot_energy_vs_rmsd(
        set1_rmsds_x, set1_energies_y, set2_rmsds_x, set2_energies_y,
        x_ticks, x_label, y_label, point_size, set1_color, set2_color,
        set1_marker, set2_marker, output_filename):
    """Plots energies for Ca-lRMSDs to the native structure on two sets
    of decoys. Saves the plot on .png and .pdf file formats.

    Args:
        set1_rmsds_x: A list/numpy array containing all the Ca-RMSD
            values for the first set.
        set1_energies_y: A list/numpy array containing all the energy
            values for the first set.
        set2_rmsds_x: A list/numpy array containing all the Ca-RMSD
            values for the second set.
        set2_energies_y: A list/numpy array containing all the energy
            values for the second set.
        x_ticks: A list/numpy array indicating ticks for X axis.
        x_label: A string indicating the label along the X axis.
        y_label: A string indicating the label along the Y axis.
        point_size: An int indicating the size of each point.
        set1_color: A string indicating the color of the points in the
            first set.
        set2_color: A string indicating the color of the points in the
            second set.
        set1_marker: A string indicating the marker of the points in the
            first set.
        set2_marker: A string indicating the marker of the points in the
            second set.
        output_filename: A string containing the name of the output
            file.
    Returns:
        None
    """

    plt.scatter(set1_rmsds_x, set1_energies_y, color=set1_color,
                s=point_size, marker=set1_marker)
    plt.scatter(set2_rmsds_x, set2_energies_y, color=set2_color,
                s=point_size, marker=set2_marker)
    plt.xlabel(x_label, fontsize=13)
    plt.ylabel(y_label, fontsize=13)
    plt.xticks(x_ticks)

    plt.savefig(output_filename + ".pdf", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)
    plt.savefig(output_filename + ".png", dpi=300, transparent='True',
                bbox_inches='tight', figsize=(5, 5), pad_inches=0.05)

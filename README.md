# Reduced Decoy Ensemble

This repository contains codes for the project detailed in the manuscript titled "Reducing Ensembles of Protein Tertiary Structures Generated De Novo via Clustering". The manuscript is currently under review and the link to the manuscript will be provided once published.

* The dataset for the codes can be found in, https://ieee-dataport.org/open-access/dataset-reducing-ensembles-protein-tertiary-structures-generated-de-novo-clustering
* The external libraries that the codes use are PyRosetta, NumPy, and Matplotlib.

## Code Files
### modules/feature_extraction.py
This module provides a function to extract the Ultrafast Shape Recognition (USR) features of a structure. It contains the following function.
* similarity_metrics: Calculates the USR features (i.e. moments of ctd, cst, fct, and ftf distances) of a given structure.
  * Arguments:
    * pose: A pyrosetta Pose object representing a structure.
  * Returns: A list containing the 12 USR features.

### modules/population.py
This module provides functions for generating initial population of structures as well as a variation operator to modify a structure. It contains the following functions.
* monte_carlo_fixed: Performs Metropolis Monte Carlo (MMC) search with a specified move, a fixed number of times in different trajectories. Each trajectory is created by performing moves from the initial structure passed in the argument.
  * Arguments:
    * pose: A pyrosetta Pose object representing initial structure
    * mover: A pyrosetta Mover object derermining the moves in MMC search.
    * score_function: A pyrosetta ScoreFunction object for scoring each move.
    * temperature: An integer/float defining the temperature of MMC search.
    * trajectory: An integer indicating the number of trajectories.
    * fixed_moves: An integer indicating the number of moves in each trajectory.
  * Returns: A list containing the population generated by the MMC search.

* monte_carlo_failure: Performs Metropolis Monte Carlo (MMC) search starting with an initial population as different trajectories. The search in each trajectory terminates when a specific number of consecutive moves fail.
  * Arguments:
    * mover: A pyrosetta Mover object derermining the moves in MMC search.
    * score_function: A pyrosetta ScoreFunction object for scoring each move.
    * temperature: An integer/float defining the temperature of MMC search.
    * population_base: A list containing the initial population.
    * successive_failures: A positive integer indicating the threshold for consecutive number of failed moves in each trajectory.
  * Returns: A list containing the population generated by the MMC search.
  
* mutation_operator: Constructs a mutation operator to perform Molecular Fragment Replacements. The operator is a pyrosetta Mover, which can be used to introduce mutation on a population or a conformation.
  * Arguments:
    * fragment_length: An integer indicating the fragment length.
    * fragment_file: A string defining the path of the file containing the fragments.
  * Returns: A pyrosetta ClassicFragmentMover object that defines each move to be a Molecular Fragment Replacement.
  
### modules/improvement.py
This module provides a function that performs local search to improve the current fitness of a given structure. It contains the following function.
* local_search: This function Performs greedy local search to improve fitness of a conformation. The local search performs specific moves to map a structure to a nearby local minimum in the energy surface. The search is terminated when a specific number of moves fail to improve the score based on a specific fitness function.
  * Arguments:
    * pose: A pyrosetta Pose object representing a structure.
    * mover: A pyrosetta Mover object that determines the moves in the local search.
    * score_function: A pyrosetta ScoreFunction object for scoring each move.
    * successive_failures: An integer indicating the threshold for consecutive number of failed moves in each trajectory.
  * Returns: A pyrosetta Pose object containing the structure with locally minimum fitness.

### modules/selection.py
This module provides a function to select next generation from current generation of individuals. It contains the following function.
* truncation: This function implements truncation selection while ensuring elitism to select a specific number of members for the next generation.
  * Arguments:
    * parent_population: A list containing members of parent population.
    * child_population: A list containing members of offspring population.
    * parents_scores: A list containing scores of each member of the parent population. The format is: [member 1 score, member 2 score, ....]. The order of members has to be consistent with parent_population argument.
    * children_scores: A list containing scores of each member of the offspring population. The format is: [member 1 score, member 2 score, ....]. The order of members has to be consistent with child_population argument.
    * elitism_rate: A float between 0 and 1 indicating the elitism percentage.
  * Returns: A list of members for the next generation of population.

### modules/io.py
This module contains the input/output utility functions needed for the project. It contains the following functions.
* get_sequence_from_fasta: This function returns the amino-acid sequence of a protein from its fasta file.
  * Arguments:
    * fasta_file: A string containing the fasta file path.
  * Returns: A string containing the amino-acid sequence of the protein.

* plot_rmsd_distribution: This function plots frequencies for Ca-lRMSDs to the native structure. Saves the plot on .png and .pdf file formats.
  * Arguments:
    * rmsds: A numpy array containing all the Ca-RMSD values.
    * x_min: A number indicating the minimum value for X axis.
    * x_max: An number indicating the maximum value for X axis.
    * x_increment: A float indicating the bin length.
    * x_label: A string indicating the label along the X axis.
    * y_label: A string indicating the label along the Y axis.
    * plot_color: A string indicating the color of the plot.
    * output_filename: A string containing the name of the output files.
  * Returns: None.

* plot_energy_vs_rmsd: This function plots energies vs. Ca-lRMSDs to the native structure of two ensembles of decoys. Saves the plot on .png and .pdf file formats.
  * Arguments:
    * set1_rmsds_x: A list/numpy array containing all the Ca-lRMSD values for the first ensemble.
    * set1_energies_y: A list/numpy array containing all the energy values for the first ensemble.
    * set2_rmsds_x: A list/numpy array containing all the Ca-lRMSD values for the second ensemble.
    * set2_energies_y: A list/numpy array containing all the energy values for the second ensemble.
    * x_ticks: A list/numpy array indicating ticks for X axis.
    * x_label: A string indicating the label along the X axis.
    * y_label: A string indicating the label along the Y axis.
    * point_size: An int indicating the size of each point.
    * set1_color: A string indicating the color of the points in the first ensemble.
    * set2_color: A string indicating the color of the points in the second ensemble.
    * set1_marker: A string indicating the marker of the points in the first ensemble.
    * set2_marker: A string indicating the marker of the points in the second ensemble.
    * output_filename: A string containing the name of the output files.
  * Returns: None.

### decoy_generation_and_feature_extraction.py
Generates decoys for a protein using HEA decoy generation algorithm and extracts Ultrafast Shape Recognition (USR) features from the decoys.

__Output:__ The output is a text file containing 14 columns and each row represents one generated decoy. The first column provides the Rosetta score4 energy, the second column provides the lRMSD to the native structure, and each of the rest of the 12 columns provides one USR feature for the structure.

__Input:__ Uses configs/generation.ini file for inputs. The inputs represented by key-value pairs in the configs/generation.ini file are explained below for each key.

* Section [init]
  * fastaPath: a string indicating the path to the .fasta file containing the amino-acid sequence of the protein.
  * pdbPath: a string indicating the path to the .pdb file containing the native tertiary structure coordinates.
  * population: an integer indicating the population size of HEA.
  * evalbudget: an integer indicating the evaluation budget of HEA.
  * noOfDecoys: an integer indicating the number of decoys to generate.

* Section [initPopulation]
  * stage1moves: an integer indicating the number of moves for the stage 1 of the initial population generation of HEA.
  * stage1score: a string indicating the name of the Rosetta energy function to use for the stage 1 of the initial population generation of HEA.
  * stage2successiveFailures: an integer indicating the number of successive failures for the stage 2 of the initial population generation of HEA.
  * stage2score: a string indicating the name of the Rosetta energy function to use for the stage 2 of the initial population generation of HEA.
  * fragmentLength: an integer indicating the length of the fragments to be used for the initial population generation of HEA.
  * fragmentFile: a string indicating the path to the file containing the fragments to use for the initial population generation of HEA.
  * stage1temperature: an integer indicating the temperature for the MMC search in the stage 1 of the initial population generation of HEA.
  * stage2temperature: an integer indicating the temperature for the MMC search in the stage 2 of the initial population generation of HEA.

* Section [improvement]
  * score: a string indicating the name of the Rosetta energy function to use for the improvement operator of HEA.
  * fragmentLength: an integer indicating the length of the fragments to be used for the improvement operator of HEA.
  * successiveFailures: an integer indicating the number of successive failures for the improvement operator of HEA.
  * fragmentFile: a string indicating the path to the file containing the fragments to use for the improvement operator of HEA.

* Section [selection]
  * elitismRate: a float between 0 and 1 indicating the elitism rate for the selection operator of HEA.
  
* Section [output]
  * filePrefix: a string indicating the path to the directory to save the output file.

### gmx_cluster.py
Applies gmx-cluster clustering on the original ensemble and generates reduced ensemble.

__Output:__ The output is a text file containing 2 columns and each row represents one structure in the reduced ensemble. The first column provides the Rosetta score4 energy and the second column provides the lRMSD to the native structure.

__Input:__ Uses configs/gmx_cluster.ini file for inputs. The inputs represented by key-value pairs in the configs/gmx_cluster.ini file are explained below for each key.

* Section [init]
  * dataInputPath: a string indicating the path to the file that contains original ensemble. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * cutoff: a float between 0 and 1 indicating the distance cutoff for the clustering.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### kmeans.py
Applies k-means clustering on the original ensemble and generates reduced ensemble.

__Output:__ The output is a text file containing 2 columns and each row represents one structure in the reduced ensemble. The first column provides the Rosetta score4 energy and the second column provides the lRMSD to the native structure.

__Input:__ Uses configs/kmeans.ini file for inputs. The inputs represented by key-value pairs in the configs/kmeans.ini file are explained below for each key.

* Section [init]
  * dataInputPath: a string indicating the path to the file that contains original ensemble. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * numClusters: an integer indicating the number of clusters for clustering.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### hierarchical.py
Applies hierarchical clustering on the original ensemble and generates reduced ensemble.

__Output:__ The output is a text file containing 2 columns and each row represents one structure in the reduced ensemble. The first column provides the Rosetta score4 energy and the second column provides the lRMSD to the native structure.

__Input:__ Uses configs/hierarchical.ini file for inputs. The inputs represented by key-value pairs in the configs/hierarchical.ini file are explained below for each key.

* Section [init]
  * dataInputPath: a string indicating the path to the file that contains original ensemble. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * numClusters: an integer indicating the number of clusters for clustering.
  * memoryPath: a string indicating the path to the directory to cache the computations.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### gmm.py
Applies GMM clustering on the original ensemble and generates reduced ensemble.

__Output:__ The output is a text file containing 2 columns and each row represents one structure in the reduced ensemble. The first column provides the Rosetta score4 energy and the second column provides the lRMSD to the native structure.

__Input:__ Uses configs/gmm.ini file for inputs. The inputs represented by key-value pairs in the configs/gmm.ini file are explained below for each key.

* Section [init]
  * dataInputPath: a string indicating the path to the file that contains original ensemble. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * numClusters: an integer indicating the number of clusters for clustering.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### truncation.py
Applies truncation on the original ensemble and generates reduced ensemble.

__Output:__ The output is a text file containing 2 columns and each row represents one structure in the reduced ensemble. The first column provides the Rosetta score4 energy and the second column provides the lRMSD to the native structure.

__Input:__ Uses configs/truncate.ini file for inputs. The inputs represented by key-value pairs in the configs/truncate.ini file are explained below for each key.

* Section [init]
  * dataInputPath: a string indicating the path to the file that contains original ensemble. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * file1path: a string indicating the path to the file that contains the reduced ensemble for hierarchical clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for k-means clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for GMM clustering.
  * file4path: a string indicating the path to the file that contains the reduced ensemble for gmx-cluster clustering.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### usr_rmsd_correlation.py
Calculates USR scores and shows correlation between the USR scores and the lRMSDs to the native structure of the decoys.

__Output:__ The outputs are Pearson's coefficient, and a plot in .png and .pdf formats.

__Input:__ Uses configs/correlation.ini file for inputs. The inputs represented by key-value pairs in the configs/correlation.ini file are explained below for each key.

* Section [init]
  * dataInputPath1: a string indicating the path to the file that contains the original ensemble for the first run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * dataInputPath2: a string indicating the path to the file that contains the original ensemble for the second run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * dataInputPath3: a string indicating the path to the file that contains the original ensemble for the third run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * dataInputPath4: a string indicating the path to the file that contains the original ensemble for the fourth run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * dataInputPath5: a string indicating the path to the file that contains the original ensemble for the fifth run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * pdbPath: a string indicating the path to the .pdb file containing the native tertiary structure coordinates.

* Section [output]
  * outputFile: a string indicating the path to the file to save the output.

### evaluation.py
Calculates and plots lRMSD distance results, plots visualization of the ensembles.

__Output:__ The outputs are size, minimum lRMSD, average lRMSD, standard deviation lRMSD of the original and the reduced ensembles. The code also outputs 4 plots in .png and .pdf formats, one for each clustering technique, to visualize the ensembles. Finally, the code outputs a plot in .png and .pdf formats to visualize the distribution of lRMSDs for the original and the reduced ensembles.

__Input:__ Uses configs/evaluate.ini file for inputs. The inputs represented by key-value pairs in the configs/evaluate.ini file are explained below for each key.

* Section [full_size]
  * file1path: a string indicating the path to the file that contains the original ensemble for the first run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * file2path: a string indicating the path to the file that contains the original ensemble for the second run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * file3path: a string indicating the path to the file that contains the original ensemble for the third run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * file4path: a string indicating the path to the file that contains the original ensemble for the fourth run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.
  * file5path: a string indicating the path to the file that contains the original ensemble for the fifth run. The file should contain 14 columns and each row should represent one generated decoy. The first column should provide the Rosetta score4 energy, the second column should provide the lRMSD to the native structure, and each of the rest of the 12 columns should provide one USR feature for the structure.

* Section [reduced_kmeans]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of k-means clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of k-means clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of k-means clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of k-means clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of k-means clustering.

* Section [reduced_hier]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of hierarchical clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of hierarchical clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of hierarchical clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of hierarchical clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of hierarchical clustering.

* Section [reduced_gmm]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of GMM clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of GMM clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of GMM clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of GMM clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of GMM clustering.


* Section [reduced_gmx]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of gmx-cluster clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of gmx-cluster clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of gmx-cluster clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of gmx-cluster clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of gmx-cluster clustering.

* Section [reduced_truncation]
  * file1path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the first run.
  * file2path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the second run.
  * file3path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the third run.
  * file4path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the fourth run.
  * file5path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the fifth run.
  
* Section [output]
  * fileName: a string indicating the path to the file to save the ensemble visualization output.
  * rmsdDist: a string indicating the path to the file to save the lRMSD distribution visualization output.

* Section [energyVsRMSD]
  * pointSize: an integer indicating the size of the points for the ensemble visualization plots.
  * set1color: a string indicating the color of the points for the original ensemble.
  * set2color: a string indicating the color of the points for the reduced ensembles.
  * set1marker: a string indicating the marker of the points for the original ensemble.
  * set2marker: a string indicating the marker of the points for the reduced ensembles.

### energy_distance_evaluation.py
Calculates and plots energy distance results.

__Output:__ The outputs are size, minimum energy distance, average energy distance, standard deviation energy distance of the original and the reduced ensembles. The code also outputs a plot in .png and .pdf formats to visualize the distribution of energy distances for the original and the reduced ensembles.

__Input:__ Uses configs/energy_distance_evaluate.ini file for inputs. The inputs represented by key-value pairs in the configs/energy_distance_evaluate.ini file are explained below for each key.

* Section [native]
  * filePath: a string indicating the path to the .pdb file containing the native tertiary structure coordinates.

* Section [full_size]
  * file1path: a string indicating the path to the file that contains the original ensemble for the first run. 
  * file2path: a string indicating the path to the file that contains the original ensemble for the second run. 
  * file3path: a string indicating the path to the file that contains the original ensemble for the third run.
  * file4path: a string indicating the path to the file that contains the original ensemble for the fourth run.
  * file5path: a string indicating the path to the file that contains the original ensemble for the fifth run. 

* Section [reduced_kmeans]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of k-means clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of k-means clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of k-means clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of k-means clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of k-means clustering.

* Section [reduced_hier]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of hierarchical clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of hierarchical clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of hierarchical clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of hierarchical clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of hierarchical clustering.

* Section [reduced_gmm]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of GMM clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of GMM clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of GMM clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of GMM clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of GMM clustering.


* Section [reduced_gmx]
  * file1path: a string indicating the path to the file that contains the reduced ensemble for the first run of gmx-cluster clustering.
  * file2path: a string indicating the path to the file that contains the reduced ensemble for the second run of gmx-cluster clustering.
  * file3path: a string indicating the path to the file that contains the reduced ensemble for the third run of gmx-cluster clustering. 
  * file4path: a string indicating the path to the file that contains the reduced ensemble for the fourth run of gmx-cluster clustering. 
  * file5path: a string indicating the path to the file that contains the reduced ensemble for the fifth run of gmx-cluster clustering.

* Section [reduced_truncation]
  * file1path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the first run.
  * file2path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the second run.
  * file3path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the third run.
  * file4path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the fourth run.
  * file5path: a string indicating the path to the file that contains the reduced ensemble generated by truncation for the fifth run.
  
* Section [output]
  * energyDist: a string indicating the path to the file to save the energy distance distribution visualization output.

### no_of_clusters.py
Plots the distributions of the number of clusters for the clustering techniques.

__Output:__ The code outputs 4 plots in .png and .pdf formats to visualize the distribution of the number of clusters, one for each clustering technique.

__Input:__ Uses configs/no_of_clusters.ini file for inputs. The inputs represented by key-value pairs in the configs/no_of_clusters.ini file are explained below for each key.

* Section [reduced_hier]
  * file1path: a string indicating the path to the file that contains the number of clusters for all the runs of hierarchical clustering on all the targets.
  * file2path: a string indicating the path to the file that contains the number of clusters for all the runs of k-means clustering on all the targets.
  * file3path: a string indicating the path to the file that contains the number of clusters for all the runs of GMM clustering on all the targets.
  * file4path: a string indicating the path to the file that contains the number of clusters for all the runs of gmx-cluster clustering on all the targets.
  
* Section [distribution]
  * color: a string indicating the color of the distribution plots.

* Section [output]
  * fileName: a string indicating the path to the file to save the number of cluster distribution output.

### bars.py
Generates bar charts for the results of the lRMSDs to the native structure for each CASP target proteins.

__Output:__ The code outputs a plot in .png and .pdf formats to visualize minimum, average, and standard deviation of lRMSDs for the original and the reduced ensembles.

__Input:__ Uses configs/bars.ini file for inputs. The inputs represented by key-value pairs in the configs/bars.ini file are explained below for each key.

* Section [input]
  * inputDir: a string indicating the directory to the files that contain the information needed to plot the bar charts. The directory contains one file for each target. Each file should contain 6 rows that provide the lRMSD value for original ensemble, reduced ensemble for hierarchical clustering, reduced ensemble for k-means clustering, reduced ensemble for GMM clustering, reduced ensemble for gmx-cluster clustering, and reduced ensemble for truncation, respectively.

* Section [output]
  * outputDir: a string indicating the path to the directory to save the bar charts.

# Developed by: Shehu Lab
# Generates decoys and extracts USR features from the decoys

from modules import population as pl
from modules import io
from modules import improvement as imp
from modules import selection as st
from modules import feature_extraction as fe
import pyrosetta as pr
import configparser as cp
import math
import random

# ***************
# Initial setup *
# ***************
pr.init()

# Read configuration file
conf = cp.ConfigParser()
conf.read("configs/generation.ini")

# Read fasta file and generate a pose from the sequence
pose = pr.pose_from_sequence(
    io.get_sequence_from_fasta(conf['init']['fastaPath'])
)

# Native pose
native_pose = pr.pose_from_pdb(conf['init']['pdbPath'])

# Convert the pose to coarse-grained representation
switch = pr.SwitchResidueTypeSetMover("centroid")
switch.apply(pose)
switch.apply(native_pose)

# Create score function for evaluation
score4_function = pr.create_score_function("score3", "score4L")
score4_function.set_weight(pr.rosetta.core.scoring.rama, 1)

population_size = int(conf['init']['population'])

eval_budget = int(conf['init']['evalbudget'])
no_of_decoys = int(conf['init']['noOfDecoys'])

min_ca_rmsd = math.inf
lowest_energy = math.inf

min_ca_rmsd_pose = pr.Pose()
lowest_energy_pose = pr.Pose()

# *****************************
# Generate initial population *
# *****************************

pop = pl.Population(native_pose)

# Create score function for stage1 MMC
scorefxn = pr.create_score_function(conf['initPopulation']['stage1score'])

# Create MFR mover for initial population
mover = pop.mutation_operator(
    int(conf['initPopulation']['fragmentLength']),
    conf['initPopulation']['fragmentFile']
)

# Run stage1 MMC
stage1_pop = pop.monte_carlo_fixed(
    pose, mover, scorefxn, int(conf['initPopulation']['stage1temperature']),
    population_size, int(conf['initPopulation']['stage1moves'])
)

# Create score function for stage2 MMC
scorefxn = pr.create_score_function(conf['initPopulation']['stage2score'])

# Run stage2 MMC
parents = pop.monte_carlo_failure(
    mover, scorefxn, float(conf['initPopulation']['stage2temperature']),
    stage1_pop, int(conf['initPopulation']['stage2successiveFailures'])
)

eval_budget -= pop.total_energy_evals

# ****************************
# Run evolutionary framework *
# ****************************

improvement = imp.Improvement(native_pose)

# Create score function for local search
scorefxn = pr.create_score_function(conf['improvement']['score'])

# Get all parameters for local search and selection
frag_length = int(conf['improvement']['fragmentLength'])
frag_file = conf['improvement']['fragmentFile']
successive_failures = int(conf['improvement']['successiveFailures'])
elitism_rate = float(conf['selection']['elitismRate'])

# Create MFR mover for local search
mover = pop.mutation_operator(frag_length, frag_file)

file_prefix = conf['output']['filePrefix']
file_suffix = str(random.randint(1, 10000))
print("File Index:" + file_suffix)

f = open(file_prefix + "data/" + file_suffix + "data.txt", "a")

decoys_generated = 0
terminate = False

while eval_budget > 0:
    # Generate children and perform local search
    children = []
    parents_scores = []
    children_scores = []

    for parent in parents:
        child = pr.Pose()
        child.assign(parent)

        # Perform mutation
        mover.apply(child)

        # Perform local search
        improved_child = pr.Pose()
        improved_child.assign(improvement.local_search(
            child, mover, scorefxn, successive_failures
        ))
        children.append(improved_child)
        eval_budget -= improvement.last_op_energy_evals

        # Evaluation bookkeeping for parents
        score4_energy = score4_function(parent)
        parents_scores.append(score4_energy)

        if score4_energy < lowest_energy:
            lowest_energy = score4_energy
            lowest_energy_pose.assign(parent)

        # Calculate metrics and populate decoy pool for parents
        if score4_energy < 0:
            pose_ca_rmsd = pr.rosetta.core.scoring.CA_rmsd(native_pose, parent)
            if pose_ca_rmsd < min_ca_rmsd:
                min_ca_rmsd = pose_ca_rmsd
                min_ca_rmsd_pose.assign(parent)

            sim_full_metrics = fe.similarity_metrics(parent)
            decoys_generated += 1

            f.write(
                "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %
                (
                    score4_energy, pose_ca_rmsd, sim_full_metrics[0],
                    sim_full_metrics[1], sim_full_metrics[2],
                    sim_full_metrics[3], sim_full_metrics[4],
                    sim_full_metrics[5], sim_full_metrics[6],
                    sim_full_metrics[7], sim_full_metrics[8],
                    sim_full_metrics[9], sim_full_metrics[10],
                    sim_full_metrics[11]
                )
            )
            if decoys_generated == no_of_decoys:
                terminate = True
                break

        # Evaluation bookkeeping for children
        score4_energy = score4_function(improved_child)
        children_scores.append(score4_energy)

        if score4_energy < lowest_energy:
            lowest_energy = score4_energy
            lowest_energy_pose.assign(improved_child)

        # Calculate metrics and populate decoy pool for children
        if score4_energy < 0:
            pose_ca_rmsd = pr.rosetta.core.scoring.CA_rmsd(native_pose,
                                                           improved_child)
            if pose_ca_rmsd < min_ca_rmsd:
                min_ca_rmsd = pose_ca_rmsd
                min_ca_rmsd_pose.assign(improved_child)

            sim_full_metrics = fe.similarity_metrics(improved_child)
            decoys_generated += 1

            f.write(
                "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %
                (
                    score4_energy, pose_ca_rmsd, sim_full_metrics[0],
                    sim_full_metrics[1], sim_full_metrics[2],
                    sim_full_metrics[3], sim_full_metrics[4],
                    sim_full_metrics[5], sim_full_metrics[6],
                    sim_full_metrics[7], sim_full_metrics[8],
                    sim_full_metrics[9], sim_full_metrics[10],
                    sim_full_metrics[11]
                )
            )

            if decoys_generated == no_of_decoys:
                terminate = True
                break

    if terminate:
        break

    # Perform selection
    parents = st.truncation(
        parents, children, parents_scores, children_scores, elitism_rate
    )

f.close()

switch = pr.SwitchResidueTypeSetMover("fa_standard")
switch.apply(lowest_energy_pose)
switch.apply(min_ca_rmsd_pose)

lowest_energy_pose.dump_pdb(file_prefix + "poses/" + file_suffix +
                            "energy.pdb")
min_ca_rmsd_pose.dump_pdb(file_prefix + "poses/" + file_suffix + "rmsd.pdb")

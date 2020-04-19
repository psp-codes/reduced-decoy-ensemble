# improvement.py
# since: 06/2018
# Developed by: Shehu Lab

"""Module for improving the fitness of a structure.

This module provides functionalities like local search to improve the
current fitness of a given structure.

Available Classes:
- Improvement: Encapsulates the operations to improve fitness of a
    conformation.
"""

import pyrosetta as pr
import math


class Improvement:
    """Encapsulates the operations to improve fitness of a conformation.

    Provides functionalities like local search to improve the current
    fitness of a given conformation.

    Public Attributes:
    - native_pose: Contains the native pose provided in the constructor
        (pyrosetta Pose object).
    - total_energy_evals: Total number of energy evaluations done in
        all the operations performed (integer).
    - last_op_energy_evals: Number of energy evaluations done in the
        last operation performed (integer).
    - min_ca_rmsd: Minimum Ca-RMSD value to the native conformation
        among all the conformations generated in all the operations
        (float)
    - last_op_min_ca_rmsd: Minimum Ca-RMSD value to the native
        conformation among all the conformations generated in the last
        operation performed (float).
    - min_ca_rmsd_pose: Conformation with minimum Ca-RMSD value to the
        native conformation among all the conformations generated in all
        the operations (pyrosetta Pose object).
    - last_op_min_ca_rmsd_pose: Conformation with minimum Ca-RMSD value
        to the native conformation among all the conformations generated
        in the last operation performed (pyrosetta Pose object).

    Available functions:
    - local_search: Performs greedy local search to improve fitness of a
        conformation.
    """

    def __init__(self, native_pose):
        """Constructor

        Args:
            native_pose: A pyrosetta Pose object containing the native
                conformation. This is used for minimum Ca-RMSD
                calculation. If you don't need this calculation, or
                don't have the native conformation, just provide a
                random Pose object.
        """
        self.native_pose = pr.Pose()
        self.native_pose.assign(native_pose)
        self.total_energy_evals = 0
        self.last_op_energy_evals = 0
        self.min_ca_rmsd = math.inf
        self.last_op_min_ca_rmsd = math.inf
        self.min_ca_rmsd_pose = pr.Pose()
        self.last_op_min_ca_rmsd_pose = pr.Pose()

    def local_search(self, pose, mover, score_function, successive_failures):
        """Performs greedy local search to improve fitness of a
        conformation.

        This local search performs specific moves to map a conformation
        to a nearby local minimum in the energy surface. The search is
        terminated when a specific number of moves fail to improve the
        score based on a specific fitness function.

        Args:
            pose: A pyrosetta Pose object containing initial
                conformation.
            mover: A pyrosetta Mover object determining the moves in
                local search.
            score_function: A pyrosetta ScoreFunction object for scoring
                each move.
            successive_failures: An int indicating the threshold for
                consecutive number of failed moves in each trajectory.

        Returns:
            A pyrosetta Pose object containing the conformation with
            locally minimum fitness.
        """
        local_minima = pr.Pose()
        local_minima.assign(pose)

        new_pose = pr.Pose()
        new_pose.assign(pose)

        self.last_op_min_ca_rmsd = pr.rosetta.core.scoring.CA_rmsd(
            self.native_pose, new_pose
        )

        local_minima_score = score_function(local_minima)
        self.last_op_energy_evals = 1

        failed = 0
        # Perform greedy local search
        while failed < successive_failures:
            mover.apply(new_pose)

            pose_ca_rmsd = pr.rosetta.core.scoring.CA_rmsd(
                self.native_pose, new_pose
            )
            if pose_ca_rmsd < self.last_op_min_ca_rmsd:
                self.last_op_min_ca_rmsd = pose_ca_rmsd
                self.last_op_min_ca_rmsd_pose.assign(new_pose)

            current_score = score_function(new_pose)
            self.last_op_energy_evals += 1

            if current_score < local_minima_score:
                local_minima.assign(new_pose)
                local_minima_score = current_score
                failed = 0
            else:
                failed += 1

        # Bookkeeping
        self.total_energy_evals += self.last_op_energy_evals
        if self.last_op_min_ca_rmsd < self.min_ca_rmsd:
            self.min_ca_rmsd = self.last_op_min_ca_rmsd
            self.min_ca_rmsd_pose.assign(self.last_op_min_ca_rmsd_pose)

        return local_minima

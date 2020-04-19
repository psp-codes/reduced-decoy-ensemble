# selection.py
# since: 10/2018
# Developed by: Shehu Lab

"""Module for selecting next generation from current generation.

This module provides methods to select next generation from
current generation.

Available Functions:
- truncation: Selects next generation via elitism truncation selection.
"""


def truncation(parent_population, child_population, parents_scores,
               children_scores, elitism_rate):
    """Selects next generation using elitism truncation selection.

    This function implements truncation selection while ensuring elitism
    to select a specific number of members for the next generation.

    Args:
        parent_population: A list containing members of parent
            population.
        child_population: A list containing members of offspring
            population.
        parents_scores: A list containing scores of each member of the
            parent population. The format is:
                [member 1 score, member 2 score, ....]
            The order of members has to be consistent with
            parent_population argument.
        children_scores: A list containing scores of each member of the
            offspring population. The format is:
                [member 1 score, member 2 score, ....]
            The order of members has to be consistent with
            child_population argument.
        elitism_rate: A float indicating the elitism percentage.

    Returns:
        A list of members for the next generation of population.
    """
    population_size = len(parent_population)
    population_indices = list(range(population_size))

    sorted_parents_indices = [x for _, x in sorted(zip(
        parents_scores, population_indices
    ))]
    sorted_parents_scores = sorted(parents_scores)

    # Slice parent population using elitism rate
    slice_index = int(population_size * elitism_rate)
    selected_parents_indices = sorted_parents_indices[:slice_index]
    selected_parents = [parent_population[i] for i in selected_parents_indices]

    combined_population = selected_parents + child_population
    combined_scores = sorted_parents_scores[:slice_index] + children_scores
    combined_population_indices = list(range(len(combined_population)))

    sorted_population_indices = [x for _, x in sorted(zip(
        combined_scores, combined_population_indices
    ))]

    selected_population_indices = sorted_population_indices[:population_size]

    # Truncate and return
    return [combined_population[i] for i in selected_population_indices]

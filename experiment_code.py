#!/usr/bin/python3

"""
This is the code that was run on a research computer of the School of Mathematics and Statistics at The University of Sydney.
"""

import pandas as pd
import time
import sys
from triangulations_optimised import *

trials = sys.argv[1]
n = int(sys.argv[2])
times = [time.time()]

disconnected = {}
distribution = {}
avg_symmetries = {}
avg_var_degree_sequence = {}
records = pd.DataFrame(columns=["sample", "n", "genus", "symmetries", "colour_swaps", "var_degree_sequence", "weight"])

n = int(n)
disconnected[n] = 0
distribution[n] = [0 for g in range(int((n-1)/2) + 1)]
avg_symmetries[n] = [0 for g in range(int((n-1)/2) + 1)]
avg_var_degree_sequence[n] = [0 for g in range(int((n-1)/2) + 1)]

tt_genus = 0
tt_gem = 0
tt_symmetry = 0
tt_weight = 0
searching = 0
sampling = 0

ttt = time.time()
for i in range(int(trials)):
    
    tt = time.time()
    while True:
        tttt = time.time()
        mu = partition_to_permutation(sample_partition(n))
        sigma = sample_permutation(n)
        sampling += time.time() - tttt

        # here: check connectedness and number of cycles...


        T = Triangulation(mu, sigma)

        if T.connected_graph():
            break
        else:
            # Disconnected graphs
            disconnected[n] += 1
    searching += time.time() - tt


    tt = time.time()
    T.compute_genus()
    tt_genus += time.time() - tt

    tt = time.time()
    T.compute_gem()
    tt_gem += time.time() - tt

    tt = time.time()
    T.compute_num_symmetries()
    tt_symmetry += time.time() - tt

    tt = time.time()
    T.compute_weight()
    tt_weight += time.time() - tt
    
    # Genus distribution
    distribution[n][T.genus] += 1/T.weight

    # Symmetries and colour swaps
    symmetries = sum(T.num_symmetries)
    avg_symmetries[n][T.genus] += 1/T.weight * symmetries
    
    colour_swaps = 0
    for num in T.num_symmetries:
        if num != 0:
            colour_swaps += 1
    
    # Variance of degree sequence
    v = T.mu.cycles + T.sigma.cycles + T.phi.cycles
    avg_degree = 6*n / v
    sum_of_squares = 0
    
    for j in T.mu.cycle_structure:
        sum_of_squares += (2*j - avg_degree)**2 * T.mu.cycle_structure[j]
    for j in T.sigma.cycle_structure:
        sum_of_squares += (2*j - avg_degree)**2 * T.sigma.cycle_structure[j]
    for j in T.phi.cycle_structure:
        sum_of_squares += (2*j - avg_degree)**2 * T.phi.cycle_structure[j]
        
    var_degree_sequence = sum_of_squares
    avg_var_degree_sequence[n][T.genus] += 1/T.weight * var_degree_sequence
    
    # Records
    #records = records.append({"sample": i,
    #                          "n": T.n,
    #                          "genus": T.genus,
    #                          "symmetries": symmetries,
    #                          "colour_swaps": colour_swaps,
    #                          "var_degree_sequence": var_degree_sequence,
    #                          "weight": T.weight}, ignore_index=True)

# Correct the averages by dividing by the sum of the inverse weights stored in distribution[n][g]
for g in range(int((n-1)/2) + 1):
    if distribution[n][g] != 0:
      avg_symmetries[n][g] = avg_symmetries[n][g] / distribution[n][g]
      avg_var_degree_sequence[n][g] = avg_var_degree_sequence[n][g] / distribution[n][g]
    else:
      avg_symmetries[n][g] = None
      avg_var_degree_sequence[n][g] = None

# Print time taken in minutes
times.append(time.time())
print("n = {}: {} minutes".format(n, round((times[-1] - times[-2]) / 60, 2)) )

# Save records and results to file
results = pd.DataFrame(data={"disconnected": disconnected,
                             "distribution": distribution,
                             "avg_symmetries": avg_symmetries, 
                             "avg_var_degree_sequence": avg_var_degree_sequence})

results.to_pickle("results_trials_{}_nlist_{}".format(trials, str(n)))
#records.to_pickle("records_trials_{}_nlist_{}".format(trials, str(n)))

everything = time.time() - ttt

print("\n")
print("sampling:             ",sampling)
print("failed samples:       ",searching)
print("genus computation:    ",tt_genus)
print("gem computation:      ",tt_gem)
print("symmetry computation: ",tt_symmetry)
print("weight computation:   ",tt_weight)
print("everything else:      ",everything-tt_genus-tt_gem-tt_symmetry-tt_weight-searching)

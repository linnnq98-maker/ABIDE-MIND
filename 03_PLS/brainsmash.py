
# -*- coding: utf-8 -*-
"""
BrainSMASH surrogate map generation
"""
import numpy as np
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.eval import sampled_fit
from brainsmash.mapgen.sampled import Sampled

# --------------------------
# File paths
# --------------------------

coord_file = r"path/to/data/L_hemisphere_regions_coordinates.txt"
output_dir = r"path/to/output/fun_brainmap"
brain_map_file = r"path/to/data/NC_group_mean_first_gradient.txt"
# --------------------------
# Load brain map
# --------------------------

brain_map = np.loadtxt(brain_map_file)
# --------------------------
# Generate spatial distance matrix
# --------------------------
filenames = volume(coord_file, output_dir)
# --------------------------
# Parameter settings
# --------------------------
N = len(brain_map)
ns_target = 490
knn_target = 1380

ns_safe = min(ns_target, N)
knn_safe = min(knn_target, N - 1)

kwargs = {
    'ns': ns_safe,
    'knn': knn_safe,
    'pv': 40,
    'nh': 6,
    'kernel': 'uniform'
}

print(
    f"Parameters: ns={ns_safe}, "
    f"knn={knn_safe}, N={N}"
)
# --------------------------
# Evaluate spatial-autocorrelation fit
# --------------------------
sampled_fit(
    brain_map,
    filenames['D'],
    filenames['index'],
    nsurr=10,
    **kwargs
)
# --------------------------
# Generate surrogate maps
# --------------------------
gen = Sampled(
    x=brain_map,
    D=filenames['D'],
    index=filenames['index'],
    **kwargs
)

surrogate_maps = gen(n=2000)
# --------------------------
# Save surrogate maps
# --------------------------
np.savetxt(
    r"path/to/output/NC_surrogate_gradient.csv",
    surrogate_maps,
    delimiter=","
)
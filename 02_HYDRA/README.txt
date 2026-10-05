# HYDRA Subtyping Analysis

This folder contains the input data, parameter settings, and documentation
used for HYDRA (Heterogeneity through Discriminative Analysis) subtyping
analysis.

## Input data

The input file used for HYDRA was:

`WND_matrix_combated_with_label.csv`

The data matrix was organized as follows:

- First column: subject ID
- Intermediate columns: MIND-derived imaging features
- Last column: group label
  - `1`: ASD
  - `-1`: control

## HYDRA analysis

HYDRA was performed in MATLAB using the following command:

```matlab
hydra('-i','WND_matrix_combated_with_label.csv','-o','.','-k',10,'-f',5);
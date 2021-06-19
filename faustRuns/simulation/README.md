# Simulation settings

To reproduce the figure in the manuscript, the files `num_of_clusters_sim.R` and `cvauc_simulation.R` must be run twice.

Run them in their current settings will produce the multivariate Gaussian output.

Then, switch the parameter `gaussianSwitch` to `FALSE` in both files. Re-run, to produce simluated data transformed to resemble CyTOF data.

These files take some time to run to completion.
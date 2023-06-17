# state_space_plot
Matlab code for dimensionality reduction and visualization of multivariate responses.
The code was designed for neural LFP recordings from multiple electrodes, but it can be used for other data-types too. Visualization is done either as ‘snapshots’ (each time-point separately) or ‘state-space trajectories’ (all time-points together).


This was written as part of the analysis for the following manuscript (so please cite if you use it): Vishne et al., Cell Reports 2023, 'Distinct Ventral Stream and Prefrontal Cortex Representational Dynamics during Sustained Conscious Visual Perception' (biorxiv DOI, to be updated when formally published): https://doi.org/10.1101/2022.08.02.502469. Example output from this repository can be seen in Figures 2A and Supplementary Figures 2A and 11A.



There are two core functions ‘do_dim_red.m’ and ‘plot_mv_patterns.m’ performing dimensionality reduction and visualization, respectively. The other functions in the folder are called from inside these functions.


Dimensionality reduction (do_dim_red.m): wrapper to Matlab’s various dimensionality reduction functions (pca\ tsne\ mds or cmds).
- The input to this function is segmented data of the form [conditions\trials x electrodes\recording-channels x time], which is then reduced on the electrodes dimension and passed to plot_mv_patterns.m for plotting.
- For non-time-resolved data simply insert data as [conditions x channels] and use the ‘snapshot’ mode.


Visualization (plot_mv_patterns.m), which has two main formats:
1. ‘Snapshots’ - each time-point is visualized in a separate plot. Dimensionality reduction can be done to all time-points together, using the previous time-point as a seed, or completely independently.
2. ‘Trajectories’ - all time-points are plotted together. In this case dimensionality reduction must be linked.


Additionally, you can:
- Add images to each data-point \ trajectory (e.g. actual images viewed by the subjects, or something representing the specific condition plotted in each case).
- Generate a video of the response evolving through time.
Both are implemented by inserting optional inputs to plot_mv_patterns.m (inputs to the plotting func).


Gal Vishne, June 2023

Twitter: @neuro_gal

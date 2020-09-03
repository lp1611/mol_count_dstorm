# mol_count_dstorm
# This repo contains code and data generating the results in the paper "Blinking Statistics and Molecular Counting in direct Stochastic Reconstruction Microscopy (dSTORM)" by Lekha Patel, David Williamson, Dylan M. Owen and Edward A.K. Cohen.

This repo contains data (.mat) and code (.m and .c) for the above paper. 
Code generating the figures are given by files: 
- FIGURE_2_sum_of_locs_plots.m for Figure 2 
- FIGURE_3_hist_sims.m for Figure 3
- FIGURE_4_alexa_plot_all.m for Figure 4
- FIGURE_6_AF647_train_test.m for Figure 6
- FIGURE_S2_alexa_plots_pshmm_mixture.m for Figure S2 (Supplement)

Data: 
- p_sum_of_locs_sims.mat and p_sum_of_locs_true.mat; simulated and theoretical probabilities generated for Figure 2
- sims_xstate_y.mat for x=3,4,5 and y=fast,medium,slow (9 datasets) for Figure 3 
- testing_number_of_locs_alexa_2.mat (Alexa Fluor 647 localization data applied to the novel method addressed in the paper) for Figure 4 
- testing_number_of_locs_alexa_negbin.mat (Alexa Fluor 647 localization data applied to the exisiting mixed Poisson-geometric method) for Figure S2 (supplement) 
- AF647_2pc405_i_by_frame.csv for i=1:22, training dataset for the T-cell study 
- AF647_est_params.mat, estimated parameter set using the PSHMM method on the training set (above) for the T-cell study 
- AF647DriCor_i.csv for i=1:4, Drift corrected T-cell testing data 
- coords.txt, central coordinates for 3x3 nano metre square regions of the cell where the our method was applied 
- testing_number_of_locs_AF647.mat, results of our method applied to each T-cell image above





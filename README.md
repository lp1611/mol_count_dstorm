# mol_count_dstorm
# This repo contains code and data generating the results in the paper "Blinking Statistics and Molecular Counting in direct Stochastic Reconstruction Microscopy (dSTORM)" by Lekha Patel, David Williamson, Dylan M. Owen and Edward A.K. Cohen.

This repo contains data (.mat) and code (.m and .c) for the above paper. 
Code generating the figures are given by files: 
- FIGURE_1_sum_of_locs_plots.m for Figure 1 
- FIGURE_3_hist_sims.m for Figure 2
- FIGURE_3_alexa_plot_all.m for Figure 3
- FIGURE_4_AF647_train_test.m for Figure 4
- FIGURE_S1_alexa_dens_plots_pshmm.m for Figure S1 (Supplement)
- FIGURE_S2_alexa_plots_pshmm_all_mixture.m for Figure S2 (Supplement)

Algorithms:
- Algorithm 2 is given by emission_delta_endstate_istate.m where i=3,4,5 for the 3,4,5 state cases (m=0,1,2 respectively). 
- Algorithm 3 is given by no_locs_pmf.m 
- Algorithm 4 is provided in FIGURE_6_AF647_train_test.m (in the final for loop).

Data: 
- p_sum_of_locs_sims.mat and p_sum_of_locs_true.mat; simulated and theoretical probabilities generated for Figure 1
- sims_xstate_y.mat for x=3,4,5 and y=fast,medium,slow (9 datasets) for Figure 2
- testing_number_of_locs_alexa_2.mat (Alexa Fluor 647 localization data applied to the novel method (PSHMM) addressed in the main text  (d=2)) for Figure 3 
- testing_number_of_locs_alexa_3state.mat (Alexa Fluor 647 localization data applied to PSHMM model with d=0) for Figure S2 (supplement) 
- testing_number_of_locs_alexa_4state.mat (Alexa Fluor 647 localization data applied to PSHMM model with d=1) for Figure S2 (supplement) 
- testing_number_of_locs_alexa_negbin.mat (Alexa Fluor 647 localization data applied to the exisiting mixed Poisson-geometric method) for Figure S2 (supplement) 
- AF647_2pc405_i_by_frame.csv for i=1:22, training dataset for the T-cell study (Figure 4)
- AF647_est_params.mat, estimated parameter set using the PSHMM method on the training set (above) for the T-cell study (Figure 4)
- AF647DriCor_i.csv for i=1:4, Drift corrected T-cell testing data (Figure 4)
- coords.txt, central coordinates for 3x3 nano metre square regions of the cell where the our method was applied (Figure 4)
- testing_number_of_locs_AF647.mat, results of our method applied to each T-cell image above (Figure 4)





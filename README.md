# extreme_value_theory
Tools of extreme value theory for analysing tail behavior of data sets.
Note, unless otherwise specified, all calculations use input files that contain the raw data in log10 and in a .pkl-file. The pickling is done by executing flux_log_quants_to_file.py.

File contents:
-raw data in ViRBO_ULF.txt and eflux.txt
-calc_errs.py contains scripts for declustering data and autocorrelation of a time series aobe a threshold value u
-declustering_autoc.py declusters the raw data and calculates the autocorrelation of the delustered maxima above a threshold u
-flux_quants_to_file.py does the same as flux_log_quants_to_file.py, but converts the data values first to measurement values by raising them to a power of 10
-fat_tail_plots.py does qq, maxiomumn to sum, mean excess and zipf plots of the raw data (doesn't use pickled data)
-fat_tail_stats.py contains function for e.g. ordering the data, fitting a GPD, calculating PIckand and Hill estimators etc. (these are used by fat_tail_plots.py)
-main_MI.py contains an object class with methods e.g. for calculating mutual information, and some of the methods are used in fat_tail_stats.py for manipulating data
-aux.py contains auxiliary functions for main_MI.py
-plot_GPD_to_log_data_bb_new.py fits a GPD to the data after declustering it. The function doens't use .pkl data. Executing the file does the fitting and saves the produced figures. The file also contains a script for doing the declusgtering using Bayesian blocks, but combining the tools was not succesful here.
-BayesBlocks.py contains the algorithm for the Bayesin block methodology (Jeffrey Scargle et al.)

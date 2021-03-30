p-value combination for PCAWG drivers

Note: we are using parts of the original code for the Empirical Brown's method (Poole et al. Bioinformatics, 2016). For details, check the following github repo: https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM

Requirements: R-3.4, uses package reshape.

The inputs to the R script are: tissue, target, dir, sif_filepath. The outputs contain a a diagnostic QQ plot of observed vs. expected p-values for all input p-value sets, a table with combined p-values, a report that indicates which p-value sets were included and a table of combined p-values after removing some outlier methods (.combined_p_values.automatic_method_removal.txt). For details please see Rheinbay, Nielsen, Abascal, Wala, Shapira et al.

aaaaaaaaa

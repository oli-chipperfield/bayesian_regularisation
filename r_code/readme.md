# Bayesian linear models, Bayesian ridge and Bayesian LASSO:
## A guide to the R code.

To make maintaining and adapting the code easier, the analysis was split across multiple files.  This document gives a brief summary of the file structure.

The analysis for each data set involves a chain of R files:

`data_prep.R`
\\
`execute_models.R`
\\
`execute_kfold_cross_validation.R`
\\
`execute_prior_tests.R`

Which executes the code for building the linear models.  The latter three files source from the files: 

`blm_stan_functions.R`
\\
`ridge_stan_functions.R`
\\
`lasso_stan_functions.R`

These files build the R functions and generated the `.stan` files, `blm.stan`, `ridge.stan` and `lasso.stan` respectively.  The execution files than saved RDS objects containing the stan model information (not included in the upload).

The analysis of the resulting model outputs was also split across a chain of R files:

`evaluate_model_parameters.R`
\\
`evaluate_cross_validation.R`
\\
`evaluate_convergence.R`

The project document makes reference to four possible data sets.  These sets of data were analysed independently in separate chains of files where the only difference was in the file roots they sourced and the name of the RDS objects they saved. The full chain for two of these data sets are defined in the R file:

`run_all_scripts.R`

The different data sets were denoted by additional underscores to the file names.  _basic file names denoted the analysis using only the original 13 parameters.  File names without appended descriptions denote the 22 variable data set where there were squared terms for all continuous variables.
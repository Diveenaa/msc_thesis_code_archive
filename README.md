# Public Health Impact of Tropical Cyclones in India using Bayesian Inference

This repositary contains the code for the MSc Computing Individual Project.
- Author: Diveena Nanthakumaran
- Supervisors: César Quilodrán Casas, Robbie M. Parks

## Dependencies

1. Causal model using Markov Chain Monte Carlo (MCMC) Sampling:

The code is run in `python` and requires the `pystan` library which can be downloaded from:
https://pystan2.readthedocs.io/en/latest/getting_started.html 

2. Bayesian ensemble model using Integrated Nested Laplace Approximation (INLA):

This code in run in `R` and requires the `R-INLA` package:
```
install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'), dependencies = TRUE)
```
```
install.packages(c('tidyverse', 'glue', 'logging', 'matrixStats'))
```

## Credits
This code has been adapted for use based on:

1. Nethery RC, Katz-Christy N, Kioumourtzoglou MA, Parks RM, Schumacher A, Anderson GB. Integrated causal-predictive machine learning models for tropical cyclone epidemiology. Biostatistics. 2023;24(2):449-464. doi:10.1093/biostatistics/kxab047 
2. Kontis V, Bennett JE, Parks RM et al. Lessons learned and lessons missed: impact of the coronavirus disease 2019 (COVID-19) pandemic on all-cause mortality in 40 industrialised countries and US states prior to mass vaccination [version 2; peer review: 2 approved]. Wellcome Open Res 2022, 6:279 (https://doi.org/10.12688/wellcomeopenres.17253.2)
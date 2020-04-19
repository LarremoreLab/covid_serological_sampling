# Estimating SARS-CoV-2 seroprevalence and epidemiological parameters with uncertainty from serological surveys

This repository stores code and data associated with the manuscript **Estimating SARS-CoV-2 seroprevalence and epidemiological parameters with uncertainty from serological surveys** by Daniel B. Larremore, Bailey K. Fosdick, Kate M. Bubar, Sam Zhang, Stephen M. Kissler, C. Jessica E. Metcalf, Caroline O. Buckee, and Yonatan H. Grad. 

This repository is maintained by [Daniel Larremore](https://larremorelab.github.io/). Questions can be directed to `daniel.larremore@colorado.edu`.

Code is written primarily in `python` but the MCMC engine has been rewritten in `R` because it is, for some reason, faster. 

# Code Availability Status

Code for this project will be made available as fast as I can document it, write docstrings, and refactor for clarity. Please pardon any delays. [4/17/2020]

# Inference of prevalence from a serological survey

When sensitivity and specificity are known, one can use those values, along with the number of positive and negative test results, to produce posterior estimates of prevalence like this one:
[![Image of Seroprevalence Posterior](https://github.com/LarremoreLab/covid_serological_sampling/blob/master/images/calculator.png){:height="50%" width="50%"}](https://larremorelab.github.io/covid-serology)

This figure can be created and downloaded by using the web-based calculator available at [https://larremorelab.github.io/covid-serology](https://larremorelab.github.io/covid-serology). 

# SEIR Simulations

In the paper, we have a figure like this one:
![Image of SEIR Simulation](https://github.com/LarremoreLab/covid_serological_sampling/blob/master/images/SEIR_workbook.png)

This figure can be recreated from scratch using the Jupyter notebook `codebase/SEIR_workbook.ipynb`. In order, this notebook:
1. Simulates hypothetical data. (Real data could be inserted here.)
2. Infers the posterior distribution over seroprevalence using the sensitivity and specificity values chosen.
3. Runs an SEIR simulation forward. (Initial conditions, parameters, and assumptions about immunity could be adjusted here.)
4. Calculates both peak timing and height distributions. (only Height is plotted, but this can be adjusted in plots.)
# Estimating SARS-CoV-2 seroprevalence and epidemiological parameters with uncertainty from serological surveys

This repository stores code and data associated with the manuscript **Estimating SARS-CoV-2 seroprevalence and epidemiological parameters with uncertainty from serological surveys** by Daniel B. Larremore, Bailey K. Fosdick, Kate M. Bubar, Sam Zhang, Stephen M. Kissler, C. Jessica E. Metcalf, Caroline O. Buckee, and Yonatan H. Grad. 

This repository is maintained by [Daniel Larremore](https://larremorelab.github.io/). Questions can be directed to `daniel.larremore@colorado.edu`.

Code is written primarily in `python` but the MCMC engine has been rewritten in `R` because it is, for some reason, faster. 

# Code Availability Status

Code for this project will be made available as fast as I can document it, write docstrings, and refactor for clarity. Please pardon any delays. [4/17/2020]

# SEIR Simulations

In the paper, we have a figure like this one:
![Image of SEIR Simulation](https://github.com/LarremoreLab/covid_serological_sampling/blob/master/images/SEIR_workbook.png)

This figure can be recreated from scratch using the Jupyter notebook `codebase/SEIR_workbook.ipynb`.
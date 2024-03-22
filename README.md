# KETCHUP
Kinetic Estimation Tool Capturing Heterogeneous Datasets Using Pyomo (KETCHUP)

# Description
Repository for Kinetic Estimation Tool Capturing Heterogeneous Datasets Using Pyomo (KETCHUP), a flexible parameter estimation tool that leverages a primal-dual interior-point algorithm to solve a nonlinear programming (NLP) problem that identifies a set of parameters capable of recapitulating the steady-state fluxes and concentrations in wild-type and perturbed metabolic networks. [1]

KETCHUP can use K-FIT [2] input files. Example K-FIT input files are located in the K-FIT repository at https://github.com/maranasgroup/K-FIT.

# Installation
In order to install directly using the environment file provided in the installation directory, anaconda needs to be previously installed. To install a self-contained environment for KETCHUP, activate the anaconda prompt and type:
```conda env create -f <path to environment.yml file>```

# Examples and Graphical User Interface
Command line python examples and a streamlit-based graphical user interface are provided. Please see https://streamlit.io for instructions on installing streamlit.

# Funding
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number DE-SC0018420). This material is based upon work supported by the Center for Bioenergy Innovation (CBI), U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number ERKP886. Funding also provided by the DOE Office of Science, Biological and Environmental Research Program Award Number DE-SC0018260. Any opinions, findings, and conclusions or recommendations expressed herein are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Computations for this research were performed on the Pennsylvania State University’s Institute for Computational and Data Sciences’ Roar supercomputer.

# References
<ol>
 <li> Mengqi Hu, Patrick F. Suthers, Costas D. Maranas. KETCHUP: parameterizing of large-scale kinetic models using multiple datasets with different reference states. Metab Eng, 2024, https://doi.org/10.1016/j.ymben.2024.02.002</li>
 <li>Saratram Gopalakrishnan, Satyakam Dash, and Costas D. Maranas, “K-FIT: An accelerated kinetic parameterization algorithm using steady-state fluxomic data,” Metab Eng, vol. 61, no. January, pp. 197–205, 2020, https://doi.org/10.1016/j.ymben.2020.03.001</li>
</ol>

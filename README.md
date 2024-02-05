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
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research under Award Number DE-SC0018420). Funding provided by The Center for Bioenergy Innovation (CBI), which is a U.S. Department of Energy Bioenergy Research Center supported by the Office of Biological and Environmental Research in the DOE Office of Science. Oak Ridge National Laboratory is managed by UT-Battelle, LLC for the US DOE under Contract Number DE-AC05-00OR22725. Any opinions, findings, and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Funding also provided by the DOE Office of Science, Office of Biological and Environmental Research (Award Number DE-SC0018260). Computations for this research were performed on the Pennsylvania State University’s Institute for Computational and Data Sciences’ Roar supercomputer.

# References
<ol>
 <li> Mengqi Hu, Patrick F. Suthers, Costas D. Maranas. KETCHUP: parameterizing of large-scale kinetic models using multiple datasets with different reference states. In preparation</li>
 <li>S. Gopalakrishnan, S. Dash, and C.D. Maranas, “K-FIT: An accelerated kinetic parameterization algorithm using steady-state fluxomic data,” Metab Eng, vol. 61, no. January, pp. 197–205, 2020, https://doi.org/10.1016/j.ymben.2020.03.001</li>
</ol>

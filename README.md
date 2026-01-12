# KETCHUP
Kinetic Estimation Tool Capturing Heterogeneous Datasets Using Pyomo (KETCHUP)

## Description
Repository for Kinetic Estimation Tool Capturing Heterogeneous Datasets Using Pyomo (KETCHUP), a flexible parameter estimation tool that leverages a primal-dual interior-point algorithm to solve a nonlinear programming (NLP) problem that identifies a set of parameters capable of recapitulating the steady-state fluxes and concentrations in wild-type and perturbed metabolic networks. [\[1\]](#cite1)

KETCHUP can use K-FIT [\[2\]](#cite2) input files. Example K-FIT input files are located in the K-FIT repository at https://github.com/maranasgroup/K-FIT.

KETCHUP is extended to time-series data [\[3\]](#cite3).

## Installation
You can download the KETCHUP repository and run it with python. 
We recommend installing it within a self-contained virtual environment using Anaconda. You can set one up from the [YAML environment file](installation/ketchup_environment.yml) in the [``installation``](installation) subdirectory by importing it using Anaconda Navigator or via the command line using:

````console
conda env create -f <path to ketchup_environment.yml file>
````

For further information, please follow the detailed installation documentation located in the [``doc`` subdirectory](doc).

## Documentation
Documentation files are located in the [``doc`` subdirectory](doc). These include a quickstart and details about the dynamic extension.

## Examples and Graphical User Interface
Command line python examples and a streamlit-based graphical user interface are provided in the [``example`` subdirectory](KETCHUP_main/example) and [``streamlit`` subdirectory](KETCHUP_main/example), respectively. README files in each location give details on usage. Streamlit is included in our suggested anaconda environment.  Please see https://streamlit.io for more details about streamlit.

## Manuscript Supplementary Materials
Files specifically associated with publications involving KETCHUP development are located in the [``Manuscript Supplementary Materials`` directory](Manuscript%20Supplementary%20Materials). These files are not required for KETCHUP to function.

## Funding
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Biological and Environmental Research under Award Number DE-SC0018420). This material is based upon work supported by the Center for Bioenergy Innovation (CBI), U.S. Department of Energy, Office of Science, Biological and Environmental Research under Award Number ERKP886. Funding also provided by the DOE Office of Science, Biological and Environmental Research Award Number DE-SC0018260. Any opinions, findings, and conclusions or recommendations expressed herein are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Computations for this research were performed on the Pennsylvania State University’s Institute for Computational and Data Sciences’ Roar supercomputer. The authors of this work recognize the Penn State Institute for Computational and Data Sciences (RRID:SCR_025154) for providing access to computational research infrastructure within the Roar Core Facility (RRID: SCR_026424).

## References
<ol>
 <li><a href="#cite1"></a>Mengqi Hu, Patrick F. Suthers, Costas D. Maranas, 2024. “KETCHUP: Parameterizing of large-scale kinetic models using multiple datasets with different reference states.” Metab. Eng., 82 (March), pp. 123–133. https://doi.org/10.1016/j.ymben.2024.02.002</li>
 <li><a href="#cite2"></a>Saratram Gopalakrishnan, Satyakam Dash, and Costas D. Maranas, 2020. “K-FIT: An accelerated kinetic parameterization algorithm using steady-state fluxomic data.” Metab. Eng., 61 (January), pp. 197–205. https://doi.org/10.1016/j.ymben.2020.03.001</li>
 <li><a href="#cite3"></a>Mengqi Hu, Bilal S. Jilani, Daniel G. Olson, and Costas D. Maranas, 2025. "Parameterization of cell-free systems with time-series data using KETCHUP." PLoS Comput Biol 21(11): e1013724. https://doi.org/10.1371/journal.pcbi.1013724</li>
</ol>

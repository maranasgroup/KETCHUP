.. KETCHUP documentation master file, created by
   sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Kinetic Estimation Tool Capturing Heterogeneous Datasets Using Pyomo (KETCHUP) is a flexible parameter estimation tool that leverages a primal-dual interior-point algorithm to solve a nonlinear programming (NLP) problem that identifies a set of parameters capable of recapitulating the steady-state fluxes and concentrations in wild-type and perturbed metabolic networks. [1] KETCHUP, by default, takes in K-FIT [2] input files and automatically constructs a kinetic model assuming mass action kinetics format explained by Tran et al., 2008 [3] and parameterizes the kinetic model with a reference flux dataset along with mutant strain flux datasets. Example K-FIT input files are located in the K-FIT repository at https://github.com/maranasgroup/K-FIT.


KETCHUP has options to automatically construct kinetic models in Michaelis-Menten format as detailed in Liebermeister and Klipp [4] and is extended to parameterize time-series data [5].


For installation instructions please refer to 'Installation.rst <https://github.com/maranasgroup/KETCHUP/doc/installation.rst>'


Welcome to KETCHUP's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   quickstart



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Funding
=======
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number DE-SC0018420). This material is based upon work supported by the Center for Bioenergy Innovation (CBI), U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number ERKP886. Funding also provided by the DOE Office of Science, Biological and Environmental Research Program Award Number DE-SC0018260. Any opinions, findings, and conclusions or recommendations expressed herein are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Computations for this research were performed on the Pennsylvania State University’s Institute for Computational and Data Sciences’ Roar supercomputer.

References 
==========
.. [1] Mengqi Hu, Patrick F. Suthers, Costas D. Maranas, 2024. “KETCHUP: Parameterizing of large-scale kinetic models using multiple datasets with different reference states.” Metab. Eng., 82 (March), pp. 123–133. https://doi.org/10.1016/j.ymben.2024.02.002
.. [2] Saratram Gopalakrishnan, Satyakam Dash, and Costas D. Maranas, 2020. “K-FIT: An accelerated kinetic parameterization algorithm using steady-state fluxomic data.” Metab. Eng., 61 (January), pp. 197–205. https://doi.org/10.1016/j.ymben.2020.03.001
.. [3] Linh M. Tran, Matthew L. Rizk, James C. Liao, 2008. "Ensemble modeling of metabolic networks." Biophys J., 95(12) (Dec). https://doi.org/10.1529/biophysj.108.135442
.. [4] Liebermeister W, Klipp E. Bringing metabolic networks to life: convenience rate law and thermodynamic constraints. Theor Biol Med Model. 2006 Dec 15;3:41. doi: 10.1186/1742-4682-3-41. PMID: 17173669; PMCID: PMC1781438.
.. [5] Mengqi Hu, Bilal S. Jilani, Daniel G. Olson, and Costas D. Maranas. "Parameterization of cell-free systems with time-series data using KETCHUP." Submitted.


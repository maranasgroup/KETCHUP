KETCHUP dynamic (extension)
===========================
KETCHUP is updated to allow for customized rate laws and parameterization fitting time-series data. We are currently still following the K-FIT input format for model and mechanism (minor update for custom rate laws) construction. Please note the data format is different than K-FIT's data input format. 

Selection of dynamic data parameterization
==========================================
Set user_options['data_type'] to dynamic and model_options['mechanism_type'] to custom

.. code-block:: python
   user_options = {
        ...previous setttings...
   'data_type' : 'dynamic'
        }
   user_options['mechanism_type'] = 'custom' # rate law to follow

Selection of data format
========================
Dynamic data input has been formated to accept .tsv files similar to COPASI dynamic data input. If the datafile is not a .tsv file KETCHUP can process excel-based data inputs. With the example dynamic data files format we developed strainer to Systmatically process the TRAINing expERimental datasets. To properly process excel-based datasets, we assume each sheet is a separate dataset following a similar initial condition format. In the KETCHUP_dynamic.py example, the data_strainer_header is such:

.. code-block:: python
   'data_strainer_header': {"t_0": ["bdh", "23bdo", "actn", "nad", "formate", "nadh", "co2"], "time": ["nadh"],
                                 "type": ["e", "c", "c", "c", "c", "c", "c", "c"],
                                 "status": ["i", "g", "g", "i", "i", "i", "i", "d"]}
where

- ``t_0`` = initial condition labels
- ``time`` = metabolite that is dependent on time measurements
- ``type`` = type of data where "e" is enzyme and "c" is component(e.g., metabolite)
- ``status`` = how to identify each label where "i" is independent data, "g" is ignored data, "d" is dependent data

In the example provided above:

- "bdh" at initial conditions is an enzyme that gives an independent measurement, this must correspond to "Enzyme ID" in the mechanism file.
- "23bd" and "actn" are ignored components (in the experiment they are not measured)
- "nad" , "formate", "nadh", "co2"  are independent measurements that only have initial conditions measured
- "time" : ["nadh"] are dependent measurements a component 

Explanation on custom rate law setup
====================================
Using formate dehydrogenase as an example where rate law is defined as:

... math::
    v(FDH) = \frac{KCAT[f] * [E] * [nad] * [hco2]}{(KI[nad]*KM[hco2] + KM[hco2]*[nad] + KM[nad]*[hco2] + [nad]*[hco2]+KI[nad]*KM[hco2]/KI[nadh]*[nadh] + KM[nad]/KM[nadh]*[hco2]*[nadh])}

In 'FDH_mechanism.xlsx' file, define the rate law (corresponding to reaction row) in the rate law column as a text file with the following format:

- KCAT[<dir>] where <dir> is either 'f' for forward or 'r' for reverse
- KM[<met>] where <met> is the metabolite for the Michaelis-Menten constant
- KI[<met>] where <met> is the metabolite for the inhibitor constant
- KCONS[<met>] where <met> is the metabolite for a kinetic constant. KCONS is an arbitrary kinetic constant that is not defined as the turnover, Michaelis-Menten, or inhibitor constant.
- [E] is by default always the enzyme concentration
- [<met>] is the metabolite concentration

KETCHUP's custom input function uses regex to automatically determine the equation in the stated format and recasts this custom rate law as a constraint.

Output file format
==================
The output results will store the kinetic parameters as a dictionary object in a .json format with the following:
kp (kinetic parameters): separated by the following

-  KCAT - turnover number for '<reaction>_<dir>' where <dir> is either 'f' or 'r'
-  KM - Michaelis-Menten constant for '<reaction>_<met>' where <met> is the metabolite
-  KI - Inhibitor constant for '<reaction>_<met>' where <met> is the metabolite
-  KCONS - Kinetic constant for '<reaction>_<met>' where <met> is the metabolite
- c (concentration of metabolites): separated by the training datasets
- e (concentration of enzymes and respective complexes): separated by the training datasets
- rate (flux rates): separated by the training datasets
- SSR (Sum of Squares residual calculated by objective function)
- time (time required to parameterize the model)
Quick start
===========

KETCHUP uses Pyomo to formulate the parameterization problem. KETCHUP uses the Pandas dataframe to transfer information required for construction of kinetic models. The tool currently uses the K-FIT input format but additional types of data formatting can be accepted by introducing a new input function.
For the automatic construction of kinetic models, user input regarding the location of model/data and their format are required.
An example script 'KETCHUP_example.py' for a steady-state problem is provided in the `example directory <../KETCHUP_main/example>`__.
To run the example script, run 'KETCHUP_example.py' with the following bash command.

.. code-block:: console
   python KETCHUP_example.py

By default it will use the file 'ipopt.opt' located in the same directory for the IPOPT solver options files, and will use the random seed of 0. Both of these can be changed from the command line:

.. code-block:: console
   python KETCHUP_example.py -s 307 -so ipopt_hsl.opt

The seed can be an integer or time to use the system clock.

You can use

.. code-block:: console
    $ python KETCHUP_example.py --help

for info about the command line parameters.

See the `README.md <../KETCHUP_main/examples/README.md>`__ file in the examples directory for more information on the examples.

Configuring KETCHUP_example.py
===============================
The first function 'ketchup_model_options()' sets up required information for kinetic model construction and parameterization. This information is passed through the model_options dictionary object. 

.. code-block:: python
    # primary items
    model_options = {
        'input_format' : 'kfit', # flag for kind of input files
        'directory_model' : f"{os.getcwd()}/data/", # location of model files
        'filename_model' : 'k-ecoli74_model.xlsx', 
        'filename_mechanism' : 'k-ecoli74_mechanism.xlsx',
        'directory_data' : f"{os.getcwd()}/data/", # location of data files
        'filename_data' : 'k-ecoli74_data.xlsx',
        'directory_output' : f"{os.getcwd()}", # location of output files
        'model_name' : 'k-ecoli74' # optional but here for easy edit. used in output
        }

Additional options provided below are set are their default.

.. code-block:: python
    # secondary items
    model_options['debug'] = False # change to True for additional runtime output
    model_options['flag_output_sbml'] = True # flag to output results to SBML file
    model_options['mechanism_type'] = 'elemental' # rate law to follow, other option is 'MM' or 'Michaelis-Menten'

Kinetic parameterizations require randomly initialized starting points. KETCHUP uses the numpy package to generate these random initial values.
Seed values are gathered through the command line argument or defaults to a value of 0 if argument is not provided. Random initial values are generated using a uniform distribution by default but a logarithmic distribution (set option to = 'log') can be selected.

.. code-block:: python
    # read seed for random number generator
    try:
        seedvalue = int(sys.argv[1])
    except (ValueError, IndexError):
        seedvalue = 0

    model_options['seedvalue'] = seedvalue 

    model_options['distribution'] = 'uniform' # distribution for initialization

IPOPT, the interior point optimizer used for parameterization of KETCHUP models uses an .opt files to determine which linear solver to use (for solving the systems of equations casted by Pyomo) and termination criteria. By default the file name is set to 'ipopt.opt' but can be provided by user through command line arguments.

.. code-block:: python
    # read filename for solver options
    try:
        model_options['solver_opt_fn'] = str(sys.argv[2])
    except (ValueError, IndexError):
        # TODO change override to use default if seed specific file does not exist
        model_options['solver_opt_fn'] = 'ipopt.opt' if int(seedvalue) <= 0 else f"{os.getcwd()}/options/ipopt_{seedvalue}.opt"

Build and parameterize KETCHUP model in KETCHUP_example.py
==========================================================
After default options for a KETCHUP run is set, we move forward with building the KETCHUP model.
First we import required packages that contain the definition for setting up the kinetic model parameterization problem.

.. code-block:: python
    import os
    from os.path import join
    import sys
    
    # add path to ktools if not installed
    sys.path.insert(0, f"{os.getcwd()}/../src/")
    
    import ktools
    from ktools.ketchup import ketchup_generate_model, solve_ketchup_model
    from ktools.io import result_dump, create_sbml_kinetic_model
    from ktools.ketchup.analysis import evaluate_stability
    from timeit import default_timer as timer

Generate the options file from previously defined function and use the options to create the kinetic model.

.. code-block:: python
    ketchup_options = ketchup_model_options()
    ketchup_model = ketchup_generate_model(ketchup_options)

Parameterization of the kinetic model with an objective function focuses on minimizing a weighted least-squares sums equation of the experimental and predicted flux rates.
The objective function is defined in ./src/ketchup/ketchup.py. Results of the parameterized model is stored in results and the status of the parameterization ('Optimal', 'infeasible', or 'MaxIterations') is stored in the status variable.

.. code-block:: python
    results = solve_ketchup_model(ketchup_model, ketchup_options)
    status = str(results.Solver[0]['Termination condition'])


The output results are stored in a dictionary object saved to a .json file. The format of the object (for elemental kinetics) is as follows:

- kf (forward kinetic parameters)
- kr (reverse kinetic parameters)
- c (concentration of metabolites): separated by the training datasets
- e (concentration of enzyme and their complexes): separated by the training datasets
- rate (flux rates): separated by the training datasets
- SSR (Sum of Squares residual calculated by objective function)
- time (time required to parameterize the model)








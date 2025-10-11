#!/usr/bin/env python3
#
# Example KETCHUP dynamic script
#
# Edit items in ketchup_user_options to select a different model, data, or other options
from argparse import Namespace


def ketchup_user_options(cmd_args: Namespace) -> dict:
    """ Creates a dictionary of user options for KETCHUP modeling that override
        default values and passes command-line arguments and program option .yaml files.

    Parameters
    ----------
    cmd_args : Namespace
        Command line arguments from argparse.

    Returns
    -------
    dict
        A dictionary containing the full set of model options.
    """
    import os
    from ktools.ketchup.ketchup import ketchup_model_options
    
    # primary items
    user_options: dict = {
        'directory_model' : os.path.join(os.getcwd(), "data"), # location of model files
        'filename_model' : 'FDH_model.xlsx',
        'filename_mechanism' : 'FDH_mechanism.xlsx',
        'mechanism_type': 'custom', # rate law to follow
        'directory_data' : os.path.join(os.getcwd(), "data"), # location of data files
        'filename_data' : 'FDH_dataset_series_A1.xlsx',
        'data_type': 'dynamic',
        'data_format': 'strainer', # Systematically processing TRAINing expERimental datasets 
        'data_strainer_header': {"t_0": ["fdh", "23bdo", "actn", "nad", "formate", "nadh", "co2"], "time": ["nadh"],
                                 "type": ["e", "c", "c", "c", "c", "c", "c", "c"],
                                 "status": ["i", "g", "g", "i", "i", "i", "i", "d"]},
        'directory_output' : os.getcwd(), # location of output files
        'model_name' : 'FDH', # optional but here for easy edit. used in output
        }

    # secondary items
    user_options['debug'] = False # default for all functions is False. change to True for additional runtime output
    #user_options['mechanism_type'] = 'custom' # rate law to follow
    user_options['distribution'] = 'uniform' # distribution for initialization
    user_options['time_delay'] = {"A1": 0, "A2": 0, "A3" : 0,
                                  "A4": 0, "A5": 0, "A6" : 0,
                                  "A7": 0, "A8": 0, "A9" : 0}
                                  # time-delay of beginning of simulations. units follow time data.
                                  # if all are the same a scalar can be also used instead of a dictionary:
                                  # user_options['time_delay'] = 0

    # options 'seedvalue' and 'filename_solver_opt' can be set/changed through command line arguments

    user_options['filename_solver_opt'] = 'ipopt.opt'

    # finally, process user defined options and command arguments to set the model options. Here,
    #   items that were not otherwise set take their default values unless overridden by
    #   the command arguments and, if set, the program-options file.

    model_options = ketchup_model_options(user_options, cmd_args)

    return model_options


def main() -> None:
    """
    Main function to set up, generate, solve, and output KETCHUP models.
    """
    import os
    import sys

    # add path to ktools if not installed
    sys.path.insert(0, os.path.join(os.getcwd(), "..", "src"))

    import ktools
    from ktools.ketchup import (ketchup_generate_model, solve_ketchup_model, ketchup_output_write,
                                ketchup_argument_parser)
    from ktools.io import result_dump, create_sbml_kinetic_model
    from ktools.ketchup.analysis import evaluate_stability, infeasible_constraints
    from timeit import default_timer as timer

    # parse command line arguments
    args = ketchup_argument_parser()

    # now read in the options for model and data files. Edit ketchup_model_options() for the
    #    specific problem, if desired. Note that any program options yaml file passed through
    #    command line arguments takes precedent, followed by command line options

    ketchup_options = ketchup_user_options(args)

    # create the model
    ketchup_model = ketchup_generate_model(ketchup_options)

    # solve the ketchup model
    time_start = timer();
    results = solve_ketchup_model(ketchup_model, ketchup_options)
    time_end = timer()

    # output results
    ketchup_output_write(results, ketchup_model, ketchup_options, time_start, time_end)

    return
    # end of main


if __name__ == '__main__':
    main()


# end of example

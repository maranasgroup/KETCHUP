#!/usr/bin/env python3
#
# Example KETCHUP script
#
# Edit items in ketchup_model_options to select a different model, data, or other options
    
def ketchup_model_options() -> dict:
    import os
    import sys
    
    # primary items
    model_options = {
        'directory_model' : f"{os.getcwd()}/data/", # location of model files
        'filename_model' : 'k-ecoli74_model.xlsx',
        'filename_mechanism' : 'k-ecoli74_mechanism.xlsx',
        'directory_data' : f"{os.getcwd()}/data/", # location of data files
        'filename_data' : 'k-ecoli74_data.xlsx',
        'directory_output' : f"{os.getcwd()}", # location of output files
        'model_name' : 'k-ecoli74' # optional but here for easy edit. used in output
        }

    # secondary items
    model_options['debug'] = False # default for all functions is False. change to True for additional runtime output

    model_options['mechanism_type'] = 'elemental' # rate law to follow

    # read seed for random number generator
    # TODO add ability to use string 'rand' for random seed or 'time' for time-based seed
    try:
        seedvalue = int(sys.argv[1])
    except (ValueError, IndexError):
        seedvalue = 0

    model_options['seedvalue'] = seedvalue 

    model_options['distribution'] = 'uniform' # distribution for initialization
    
    # read filename for solver options
    try:
        model_options['solver_opt_fn'] = str(sys.argv[2])
    except (ValueError, IndexError):
        # TODO change override to use default if seed specific file does not exist
        model_options['solver_opt_fn'] = 'ipopt.opt' if int(seedvalue) <= 0 else f"{os.getcwd()}/options/ipopt_{seedvalue}.opt"

    return model_options

def main() -> None:
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

    # first read in the model and data files. Edit ketchup_model_options() for the
    # specific problem.

    ketchup_options = ketchup_model_options()

    # create the model

    ketchup_model = ketchup_generate_model(ketchup_options)

    # solve the ketchup model

    time_start = timer();
    try:
        results = solve_ketchup_model(ketchup_model, ketchup_options)
        status = str(results.Solver[0]['Termination condition'])
    except ValueError:
        status = 'Status Error'
    time_end = timer()

    # output results

    seedvalue = ketchup_options['seedvalue']
    
    result_dump(f"{ketchup_options['directory_output']}/{ketchup_options['model_name']}_{status}_results",seedvalue,ketchup_model,time_end - time_start,status)

    # create SBML for the model at the solution and output to file
    kmodel_sbml = create_sbml_kinetic_model(ketchup_model)
    with open(f"{ketchup_options['directory_output']}/{ketchup_options['model_name']}_{status}_results_{seedvalue}.xml", 'w') as output_file:
        output_file.write(kmodel_sbml)

    return
    # end of main
    
if __name__ == '__main__':
    main()
    


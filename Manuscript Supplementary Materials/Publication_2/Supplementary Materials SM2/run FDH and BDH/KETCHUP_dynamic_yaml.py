#!/usr/bin/env python3
#
# KETCHUP script
#
# This file only uses info from the YAML options file.

def main() -> None:
    """
        Main function to set up, generate, solve, and output KETCHUP models.
    """
    import os
    import sys
    
    # add path to ktools if not installed
    sys.path.insert(0, f"{os.getcwd()}/../src/")
    
    import ktools
    from ktools.ketchup import (ketchup_generate_model, solve_ketchup_model,
                                ketchup_output_write,
                                ketchup_argument_parse)
    from ktools.ketchup.ketchup import ketchup_model_options
    from ktools.ketchup.analysis import evaluate_stability
    from timeit import default_timer as timer

    # parse command line arguments
    args = ketchup_argument_parse()

    # only use program options yaml file passed through command line arguments
    ketchup_options = ketchup_model_options({}, args)

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

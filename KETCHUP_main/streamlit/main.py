# Streamlit GUI for KETCHUP

def ketchup_model_options() -> dict:
    import os
    import sys
    
    # primary items
    model_options = {
        'input_format' : 'kfit', # flag for kind of input files
        'directory_model' : f"{os.getcwd()}/data/", # location of model files
        'filename_model' : 'model.xlsx',
        'filename_mechanism' : 'mechanism.xlsx',
        'directory_data' : f"{os.getcwd()}/data/", # location of data files
        'filename_data' : 'data.xlsx',
        'directory_output' : f"{os.getcwd()}", # location of output files
        'model_name' : 'kmodel' # optional but here for easy edit. used in output
        }

    # secondary items
    model_options['debug'] = False # change to True for additional runtime output
    model_options['flag_output_sbml'] = True # flag to output results to SBML file
    model_options['mechanism_type'] = 'elemental' # rate law to follow

    # set seed for random number generator
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

def run_streamlit_gui() -> None:
    import streamlit as st
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

    import numpy as np
    import tempfile
    import time

    # set ketchup options
    ketchup_options = ketchup_model_options()

    # create gui
    st.image('./images/header.jpg')

    ketchup_options['filename_model'] = None
    ketchup_options['filename_mechanism'] = None
    ketchup_options['filename_data'] = None

    model_file = st.file_uploader('Select Model file')
    mechanism_file = st.file_uploader('Select Mechanism file')
    data_file = st.file_uploader('Select Data file')
    ketchup_options['model_name'] = st.text_input('Model name', 'kmodel')
    
    @st.cache_data
    def random_seed_value() -> int:
        return int(str(int(time.time()))[-3:])

    seed_input = st.text_input('Seed value', random_seed_value() )
    ketchup_options['seedvalue'] = int(seed_input)
    if st.button('Start parameterization'):
    #check for files
        if ketchup_options['filename_model'] and ketchup_options['filename_mechanism'] and ketchup_options['filename_data']:
            # if any file are missing do not run
            st.write('Missing files')
        else:
            #create temporaray directory to store the file
            directory_model = './scratch/'
            if not os.path.isdir(directory_model): os.makedirs(directory_model)
            for file in [model_file, mechanism_file, data_file]:
                temp_path = os.path.join(directory_model, file.name)
                with open(temp_path, 'wb')  as f:
                    f.write(file.getvalue())
                f.close()
        
            ketchup_options['directory_model'] = directory_model
            ketchup_options['directory_data'] = directory_model
            ketchup_options['directory_output'] = directory_model
            ketchup_options['filename_model'] = model_file.name
            ketchup_options['filename_mechanism'] = mechanism_file.name
            ketchup_options['filename_data'] = data_file.name

            seedvalue = ketchup_options['seedvalue']
            print(f"random seed value {seedvalue}")

            # create the model

            ketchup_model = ketchup_generate_model(ketchup_options)
    
            # solve the ketchup model

            time_start = timer();
            with st.spinner('Running.....'):
                try:
                    results = solve_ketchup_model(ketchup_model, ketchup_options)
                    status = str(results.Solver[0]['Termination condition'])
                except ValueError:
                    status = 'Status Error'
            time_end = timer()

            status_msg = 'Unknown status, please check command prompt'
            if status == 'maxIterations': status_msg = 'Max Iterations reached, no solution found'
            if status == 'Status Error': status_msg = 'Unable to converge to a solution'
            if status == 'optimal' : status_msg = f"Optimal solution found - SSR : {round(float(ketchup_model.obj()),3)}"
        
            st.write(status_msg)
            
            result_dump(f"{ketchup_options['directory_output']}{ketchup_options['model_name']}_{status}_results",seedvalue,ketchup_model,time_end - time_start,status)


            if ketchup_options['flag_output_sbml']:
                kmodel_sbml = create_sbml_kinetic_model(ketchup_model)
                with open(f"{ketchup_options['model_name']}_{status}_results_{seedvalue}.xml", 'w') as output_file:
                    output_file.write(kmodel_sbml)
            pass

    return

if __name__ == '__main__':
    run_streamlit_gui()


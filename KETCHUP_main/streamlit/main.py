# Streamlit GUI for KETCHUP


def run_streamlit_gui() -> None:
    """
    Creates and runs a Streamlit GUI for KETCHUP parameterization and analysis.

    Returns
    -------
    None
    """
    import streamlit as st
    import tempfile
    import os
    from os.path import join
    import sys
    import time
    import pandas as pd
    from io import StringIO
    import shutil

    # add path to ktools if not installed
    sys.path.insert(0, os.path.join(os.getcwd(),"..","src"))

    import ktools
    from ktools.ketchup import (ketchup_generate_model, solve_ketchup_model, ketchup_output_write,
                                ketchup_argument_parser)
    from ktools.ketchup.ketchup import ketchup_model_options
    from ktools.ketchup.analysis import evaluate_stability
    from timeit import default_timer as timer
    from argparse import Namespace

    # create GUI
    st.image('./images/header.jpg')
    #base empty settings
    base_user_options = {
        'input_format' : 'kfit', 'directory_model' : f"{os.getcwd()}/data/", 'filename_model' : 'k-ecoli74_model.xlsx', 'filename_mechanism' : '',
        'directory_data' : f"{os.getcwd()}/data/",'filename_data' : '',
        'data_type' : 'static','directory_output' : f"{os.getcwd()}",'model_name' : 'k-model'}
    ketchup_options = ketchup_model_options(base_user_options,Namespace(program_options=None, seed=None, solver_options=None, time_delay=None) )
    #if YAML file is available, do not use file selections
    show_uploader = st.checkbox("YAML file?"); selected = False; 
    #default values
    model_name_opt = 'k-model';data_type_opt = 'static'; distr_opt = 'uniform'; seedvalue_opt = 0; mech_type_opt = 'elemental';filename_solver_opt = ''
    st.session_state.uploader_disabled = False
    if show_uploader:
        options_file = st.file_uploader("Choose options YAML file")  
        if options_file is not None:
            selected = True
            with tempfile.TemporaryDirectory() as tmpdir:
                with tempfile.NamedTemporaryFile(dir=tmpdir) as f:
                    f.write(options_file.getbuffer()); f.flush()
                    args = Namespace(program_options=os.path.join(tmpdir,f.name), seed=None, solver_options=None, time_delay=None) 
                    ketchup_options = ketchup_model_options({},args)
                    #grab YAML info
                    model_name_opt = ketchup_options['model_name']
                    seedvalue_opt = ketchup_options['seedvalue']
                    data_type_opt = ketchup_options['data_type']
                    distr_opt = ketchup_options['distribution']
                    mech_type_opt = ketchup_options['mechanism_type']
                    filename_solver_opt = ketchup_options['filename_solver_opt']
    else:
        selected=False

    model_file_opt = st.file_uploader('Select Model file')
    if selected and model_file_opt is None:  st.write(f"Model file (YAML): {ketchup_options['directory_model']}{ketchup_options['filename_model']}")

    mechanism_file_opt = st.file_uploader('Select Mechanism file')
    if selected and mechanism_file_opt is None:  st.write(f"Mechanism file (YAML): {ketchup_options['directory_model']}{ketchup_options['filename_mechanism']}")

    data_file_opt = st.file_uploader('Select Data file')
    if selected and data_file_opt is None:  st.write(f"Data file (YAML): {ketchup_options['directory_data']}{ketchup_options['filename_data']}")
    
    model_name_opt = st.text_input("Model name", value = model_name_opt)
    add_options = st.checkbox("Additional options")
    if add_options:
        data_type_labels = ["static","dynamic"]; ind1 = data_type_labels.index(data_type_opt)
        data_type_opt = st.selectbox("Data type",data_type_labels,index=ind1)
        seedvalue_opt = st.text_input('Seed value', seedvalue_opt )
        distr_labels = ['uniform','logarithmic'];ind2 = distr_labels.index(distr_opt)
        distr_opt = st.selectbox("Random initialization distribution style",distr_labels,index=ind2)
        mech_type_labels = ['elemental','michaelis-menten','custom']; ind3= mech_type_labels.index(mech_type_opt)
        mech_type_opt = st.selectbox("Mechanism type",mech_type_labels,index=ind3)
        filename_solver_opt = st.text_input("Solver options file", ketchup_options['filename_solver_opt'])

    if st.button('Start parameterization'):
    #check for files
        if ( model_file_opt is None or mechanism_file_opt is None or data_file_opt is None ) and not selected:
            # if any file are missing do not run
            st.write('Missing files')
        else:
            #assume relative directory paths
            ketchup_options['filename_model'] = os.path.join(ketchup_options['directory_model'],ketchup_options['filename_model'])
            ketchup_options['filename_mechanism'] = os.path.join(ketchup_options['directory_model'],ketchup_options['filename_mechanism'])
            ketchup_options['filename_data'] = os.path.join(ketchup_options['directory_data'],ketchup_options['filename_data'])  
            ketchup_options['directory_model'] = './'; ketchup_options['directory_data'] = './'

            #if file_uploader if selected then store the files in temporary location
            with tempfile.TemporaryDirectory() as tmpdir:
            #tmpdir = f'./sandbox_{(int(time.time() * 1000 ) % (2^32 - 1)) }'; os.makedirs(tmpdir)
                for id,file in {'filename_model':model_file_opt, 'filename_mechanism':mechanism_file_opt, 'filename_data':data_file_opt}.items():
                    if file is not None:
                        #if '.tsv' in file.name:
                        #    df = pd.read_csv(StringIO(file.getvalue().decode('utf-8')), sep='\t')
                        #    df.to_csv(os.path.join(tmpdir,file.name))
                        #else:
                        with open(os.path.join(tmpdir,file.name),'wb') as f:
                            f.write(file.getvalue()); f.close()
                        ketchup_options[id] = os.path.join(tmpdir,file.name)

                    ketchup_options['model_name'] = model_name_opt
                    ketchup_options['data_type'] = data_type_opt
                    ketchup_options['seedvalue'] = int(seedvalue_opt)
                    ketchup_options['distribution'] = distr_opt
                    ketchup_options['mechanism_type'] = mech_type_opt
                    ketchup_options['filename_solver_opt'] = filename_solver_opt
                ketchup_options = ketchup_model_options(ketchup_options) # update options
                ketchup_model = ketchup_generate_model(ketchup_options)
    
            # solve the ketchup model

            time_start = timer();
            with st.spinner('Running.....'):
                 results = solve_ketchup_model(ketchup_model, ketchup_options)
            time_end = timer()

            ketchup_output_write(results, ketchup_model, ketchup_options, time_start, time_end)

    return


if __name__ == '__main__':
    run_streamlit_gui()

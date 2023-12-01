import streamlit as st
import os
from os.path import join
import sys

# add path to ktools if not installed
sys.path.insert(0, f"{os.getcwd()}/../src/")

import ktools
from ktools.io import read_kfit_model_xlsx
from ktools.io import read_kfit_data_xlsx
from ktools.io import result_dump, create_sbml_kinetic_model, evaluate_stability
from ktools.core import *
import numpy as np
import tempfile
from pyomo.environ import *
from pyomo.dae import *
import time

##STREAMLIT GUI START
st.image('./images/header.jpg')
model_file = st.file_uploader('Select Model file')
mechanism_file = st.file_uploader('Select Mechanism file')
data_file = st.file_uploader('Select Data file')
model_name = st.text_input('Model name', 'kmodel')
#current defaults
mech_type = 'elemental'

#seed
@st.cache_data
def random_seed_value():
    return int(str(int(time.time()))[-3:])

seed_input = st.text_input('Seed value', random_seed_value() )

#start parameterization
model_missing = False; mechanism_missing = False; data_missing = False
if st.button('Start parameterization'):
    #check for file
    if model_file is None: model_missing = True
    if mechanism_missing is None: mechanism_missing = True
    if data_missing is None: data_missing = True
    # if files are missing do not run
    if model_missing or mechanism_missing or data_missing:
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

        #st.write(f"Files temporarily stored in {directory_model}")
        (m_model, mech_df) = read_kfit_model_xlsx(join(directory_model,model_file.name),
                                              join(directory_model,mechanism_file.name),
                                              mech_type=mech_type,debug=False)

        data_df = read_kfit_data_xlsx(join(directory_model,data_file.name))
        #pass the selected blocks of the data_dataframe to the constructor

        data_dict = create_data_dict(data_df)
        

        seedvalue = int(seed_input)
        print(f"random seed value {seedvalue}")
        ketchup_model = create_initial_model(m_model,mech_df,data_dict,seedvalue=seedvalue,distribution='uniform',mech_type=mech_type)

        experiments = []
        for i,d in enumerate(data_dict.keys()):
            ketchup_model.key = d
            ketchup_model.add_component(f"experiment{i}", Block(rule=create_sMB))
            experiments.append(eval(f"ketchup_model.experiment{i}"))
            print(f"dataset included: {d}")

        ketchup_model.add_component('block_list', Set(initialize = experiments))
        ketchup_model.reaction_rate = Constraint(ketchup_model.block_list, ketchup_model.experiment0.REACTIONS, rule=net_reaction_rate)

        #ketchup_model.experiment0.vf.pprint()
        if mech_type == 'elemental':
            def ref_conc_constraint(m,s):
                return m.experiment0.c[s] == 1.0
            
                #fix wild-type fluxes
            for r in data_dict['WT']:
                ketchup_model.experiment0.rate[r].fix(data_dict['WT'][r][0])
             
            ketchup_model.ref_conc_constraint = Constraint(ketchup_model.experiment0.SPECIES, rule=ref_conc_constraint)        
            ketchup_model.obj = Objective( sense=minimize,
                         expr = sum(  b.error for b in experiments )               
                             )
            ketchup_model.vf_rate = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_F, rule = elemental_vf)
            ketchup_model.vr_rate = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_R, rule = elemental_vr)
            ketchup_model.es_net  = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_F, rule = es_net_balance) 
            

        #load status results
        solver = SolverFactory('ipopt')
        option_name = 'ipopt.opt' 
        solver.options.option_file_name = option_name
        
        from timeit import default_timer as timer

        time_start = timer();
        with st.spinner('Running.....'):
            try: 
                results = solver.solve(ketchup_model, tee=True,symbolic_solver_labels=True);
                total_error = sum( b.error.value for exps in experiments for b in exps[:] )
                status = str(results.Solver[0]['Termination condition']) 
            except ValueError:
                status = 'Status Error'
                total_error = 'NA'
        time_end = timer()
        status_msg = 'Unknown status, please check command prompt'
        if status == 'maxIterations': status_msg = 'Max Iterations reached, no solution found'
        if status == 'Status Error': status_msg = 'Unable to converge to a solution'
        if status == 'optimal' : status_msg = f"Optimal solution found - SSR : {round(float(ketchup_model.obj()),3)}"
        
        st.write(status_msg)
        result_dump(f"{model_name}_{status}_results",seedvalue,ketchup_model,time_end - time_start,status)
        #evaluate_stability(mech_df,experiments,sys.argv[1],ketchup_model,time_end - time_start,status)
        #if status == "optimal":
        if True:
            kmodel_sbml = create_sbml_kinetic_model(ketchup_model)
            with open(f"{model_name}_{status}_results_{seedvalue}.xml", 'w') as output_file:
                output_file.write(kmodel_sbml)
        pass


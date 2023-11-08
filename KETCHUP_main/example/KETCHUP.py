import os
from os.path import join
import ktools
from ktools.io import read_kfit_model_xlsx
from ktools.io import read_kfit_data_xlsx
from ktools.io import result_dump, evaluate_stability
from ktools.core import *
import sys

directory_model = f'{os.getcwd()}/data/'
filename_model = "k-ecoli74_model.xlsx"
filename_mechanism = "k-ecoli74_mech.xlsx"
filename_data = "k-ecoli74_data.xlsx"
model_name = 'k-ecoli74'

mech_type = 'elemental'

(m_model, mech_df) = read_kfit_model_xlsx(join(directory_model,filename_model),
                                          join(directory_model,filename_mechanism),
                                          mech_type=mech_type,debug=False)

data_df = read_kfit_data_xlsx(join(directory_model,filename_data))

# pass the selected blocks of the data_dataframe to the constructor

data_dict = create_data_dict(data_df)
    
from pyomo.environ import *
from pyomo.dae import *
seedvalue = int(sys.argv[1])
dae_model = create_initial_model(m_model,mech_df,data_dict,seedvalue=seedvalue,distribution='uniform',mech_type=mech_type)

experiments = []
for i,d in enumerate(data_dict.keys()):
    dae_model.key = d
    dae_model.add_component(f'experiment{i}', Block(rule=create_sMB))
    experiments.append(eval(f'dae_model.experiment{i}'))
    print(f'dataset included: {d}')

dae_model.add_component('block_list', Set(initialize = experiments))
dae_model.reaction_rate = Constraint(dae_model.block_list, dae_model.experiment0.REACTIONS, rule=net_reaction_rate)

#dae_model.experiment0.vf.pprint()
if mech_type == 'elemental':
    def ref_conc_constraint(m,s):
        return m.experiment0.c[s] == 1.0
        
            #fix wild-type fluxes
    for r in data_dict['WT']:
        dae_model.experiment0.rate[r].fix(data_dict['WT'][r][0])
         
    dae_model.ref_conc_constraint = Constraint(dae_model.experiment0.SPECIES, rule=ref_conc_constraint)        
    dae_model.obj = Objective( sense=minimize,
                 expr = sum(  b.error for b in experiments )               
                     )
    dae_model.vf_rate = Constraint(dae_model.block_list, dae_model.ELEMENTALSTEP_F, rule = elemental_vf)
    dae_model.vr_rate = Constraint(dae_model.block_list, dae_model.ELEMENTALSTEP_R, rule = elemental_vr)
    dae_model.es_net  = Constraint(dae_model.block_list, dae_model.ELEMENTALSTEP_F, rule = es_net_balance) 
    


#load status results
solver = SolverFactory('ipopt')
option_name = 'ipopt.opt' if int(seedvalue) <= 0 else f'./options/ipopt_{seedvalue}.opt'
solver.options.option_file_name = option_name

from timeit import default_timer as timer

time_start = timer();
try:
    results = solver.solve(dae_model, tee=True,symbolic_solver_labels=True);
    total_error = sum( b.error.value for exps in experiments for b in exps[:] )
    status = str(results.Solver[0]['Termination condition']) 
except ValueError:
    status = 'Status Error'
    total_error = 'NA'
time_end = timer()
 
result_dump(f"{model_name}_{status}_results",sys.argv[1],dae_model,time_end - time_start,status)
evaluate_stability(mech_df,experiments,sys.argv[1],dae_model,time_end - time_start,status)



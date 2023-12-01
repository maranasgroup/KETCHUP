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
import sys

directory_model = f"{os.getcwd()}/data/"
filename_model = 'k-ecoli74_model.xlsx'
filename_mechanism = 'k-ecoli74_mech.xlsx'
filename_data = 'k-ecoli74_data.xlsx'
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

try:
    seedvalue = int(sys.argv[1])
except:
    seedvalue = 0

ketchup_model = create_initial_model(m_model,mech_df,data_dict,seedvalue=seedvalue,distribution='uniform',mech_type=mech_type)

try:
    ketchup_model.name = model_name
except:
    pass

experiments = []
for i,d in enumerate(data_dict.keys()):
    ketchup_model.key = d
    ketchup_model.add_component(f'experiment{i}', Block(rule=create_sMB))
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
option_name = 'ipopt.opt' if int(seedvalue) <= 0 else f"{os.getcwd()}/options/ipopt_{seedvalue}.opt"
solver.options.option_file_name = option_name

from timeit import default_timer as timer

time_start = timer();
try:
    results = solver.solve(ketchup_model, tee=True,symbolic_solver_labels=True);
    total_error = sum( b.error.value for exps in experiments for b in exps[:] )
    status = str(results.Solver[0]['Termination condition']) 
except ValueError:
    status = 'Status Error'
    total_error = 'NA'
time_end = timer()

# output results

result_dump(f"{model_name}_{status}_results",seedvalue,ketchup_model,time_end - time_start,status)

kmodel_sbml = create_sbml_kinetic_model(ketchup_model)
with open(f"{model_name}_{status}_results_{seedvalue}.xml", 'w') as output_file:
    output_file.write(kmodel_sbml)

if status == 'optimal':
    try:
        evaluate_stability(mech_df,experiments,seedvalue,ketchup_model,time_end - time_start,status)
    except:
        pass


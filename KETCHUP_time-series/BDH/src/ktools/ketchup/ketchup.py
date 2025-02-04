""" 
ketchup functions

collects basic functions that for working with KETCHUP specific
models

""" 
import pyomo.core.base
from typing import Union
#KETCHUP_DYNAMIC packages for transformation factory for collocation poitns
from pyomo.dae import *
from pyomo.environ import *

def ketchup_generate_model(ketchup_options: dict) -> pyomo.core.base.PyomoModel.ConcreteModel:
    """Returns pyomo Concrete Model for subsequent computations and analysis
    """

    import pyomo.environ
    from pyomo.environ import Constraint, Block, Objective
    import cobra.core
    import ktools
    from ktools.core import create_data_dict, create_initial_model
    from ktools.io import read_kfit_model_xlsx
    from ktools.io import read_kfit_data_xlsx
    from os.path import join

    
    m_model, mech_df = read_kfit_model_xlsx(join(ketchup_options['directory_model'],
                                                   ketchup_options['filename_model']),
                                              join(ketchup_options['directory_model'],
                                                   ketchup_options['filename_mechanism']),
                                              mech_type=ketchup_options['mechanism_type'],
                                              debug=ketchup_options['debug'])
    
    
    #data_df = read_kfit_data_xlsx(join(ketchup_options['directory_data'],
    #                                   ketchup_options['filename_data']))
    #KETCHUP_DYNAMIC modification
    #data_dict = create_data_dict(data_df)
    import openpyxl
    filename_workbook = join(ketchup_options['directory_data'],ketchup_options['filename_data'])
    workbook = openpyxl.load_workbook(filename=filename_workbook, data_only = True)
    data_dict = {}
    for s in workbook.sheetnames:
        ss = workbook[s]
        data = list(ss.values)
        enzyme_conc = data[0][1]
        bd23_conc = data[1][1]
        acetoin_conc = data[2][1]
        nad_conc = data[3][1]
        formate_conc = data[4][1]
        nadh_conc = data[5][1]
        time_data = data[8:]
        data_dict.update({s:{"enzyme":enzyme_conc , "nad": nad_conc, "formate":formate_conc,'acetoin':acetoin_conc, 'nadh':nadh_conc, '23bd':bd23_conc, "time":time_data }})
    
            
    # TODO: add code to set default values for secondary options if not in ketchup_options
    tde = ketchup_options['tde'] #time delay error

    ketchup_model = create_initial_model(m_model, mech_df, data_dict, seedvalue=ketchup_options['seedvalue'], tde = tde, distribution=ketchup_options['distribution'],mech_type=ketchup_options['mechanism_type'])
    #save ketchup_options into the model object
    ketchup_model.options = ketchup_options
    try:
        ketchup_model.name = str(ketchup_options['model_name'])
    except NameError:
        ketchup_model.name = 'unknown'

    experiments = []

    for i,d in enumerate(data_dict.keys()):
        ketchup_model.key = d
        if 'tde_data' in ketchup_options: 
            ketchup_model.tde = ketchup_options['tde_data'][d]['tde']
        ketchup_model.add_component(f'experiment{i}', Block(rule=ktools.core.create_dMB))
        print(f"dataset included: {d}")
        experiments.append(eval(f"ketchup_model.experiment{i}"))
    #ketchup_model.experiment0.pprint()   
    ketchup_model.add_component('block_list', pyomo.environ.Set(initialize = experiments))

    
    #try:
    #    basis_expt = str(ketchup_options['basis_id'])
    #except KeyError:
    #    basis_expt = 'WT' if 'WT' in data_dict else data_dict[0]


    #ketchup_model.experiment0.vf.pprint()
    if ketchup_options['mechanism_type'] == 'elemental':
        def ref_conc_constraint(m,s):
            return m.experiment0.c[s] == 1.0
        #KETCHUP_DYNAMIC no basis experiments for dynamic test
        #fix wild-type fluxes
        #for r in data_dict[basis_expt]:
        #    ketchup_model.experiment0.rate[r].fix(data_dict[basis_expt][r][0])
             
        #ketchup_model.ref_conc_constraint = Constraint(ketchup_model.experiment0.SPECIES, rule=ref_conc_constraint)     

        ketchup_model.obj = Objective( sense=pyomo.environ.minimize,
                     expr = sum(  b.error for b in experiments )               
                         )
        ketchup_model.obj.pprint()
        
        #KETCHUP_DYNAMIC renamed rule to the dynamic version
        for i,b in enumerate(ketchup_model.block_list):
            b.add_component(f'vf_rate{i}', Constraint(b.time, ketchup_model.ELEMENTALSTEP_F, rule = ktools.core.d_elemental_vf) )
            b.add_component(f'vr_rate{i}', Constraint(b.time, ketchup_model.ELEMENTALSTEP_R, rule = ktools.core.d_elemental_vr) )
            b.add_component(f'es_net{i}', Constraint(b.time,  ketchup_model.ELEMENTALSTEP_F, rule = ktools.core.d_es_net_balance) )
        ketchup_model.kr['NADH_decomp_0'].fix()
        ketchup_model.kr['NADH_decomp_1'].fix()
        ketchup_model.kr['NADH_decomp_0'].value = 1
        ketchup_model.kr['NADH_decomp_1'].value = 1
    #08/05/2024 CUSTOM RATE LAWS
    elif ketchup_options['mechanism_type'] == 'custom':

        ketchup_model.obj = Objective( sense=pyomo.environ.minimize,
                     expr = sum(  b.error for b in experiments )               
                         )
        for i,b in enumerate(ketchup_model.block_list):
            b.add_component('rate_law', Constraint(b.time, rule = ktools.core.fdh_rate) )
            b.add_component('decomp' , Constraint(b.time, rule = ktools.core.nadh_decomp) )
            
    elif ketchup_options['mechanism_type'] == 'michaelis-menten':
        
        ketchup_model.obj = Objective( sense=pyomo.environ.minimize,
                  expr = sum(  b.error for b in experiments )               
                      )          
        for i,b in enumerate(ketchup_model.block_list):
            b.add_component('decomp' , Constraint(b.time, rule = ktools.core.nadh_decomp) ) 
            b.add_component('reaction_rate', Constraint(b.time,                                                        
                                                         ketchup_model.experiment0.REACTIONS,
                                                         rule=ktools.core.d_net_reaction_rate) 
                            )         
             
    return ketchup_model


def solve_ketchup_model(ketchup_model: pyomo.core.base.PyomoModel.ConcreteModel,
                        ketchup_options: dict):
    """Returns solved pyomo Concrete Model using solver options file
    """
    #KETCHUP_DYNAMIC replace solver to use collocation
    
    solver = pyomo.environ.SolverFactory('ipopt')
    #solver.options.option_file_name = ketchup_options['solver_opt_fn']

    try:
        results = solver.solve(ketchup_model, tee=True,symbolic_solver_labels=True);
        #total_error = sum( b.error.value for exps in experiments for b in exps[:] )
        #results = solver.solve(ketchup_model, tee=True)
        status = str(results.Solver[0]['Termination condition']) 
    except ValueError:
        status = 'Status Error'

    return results




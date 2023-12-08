""" 
ketchup functions

collects basic functions that for working with KETCHUP specific
models

""" 
import pyomo.core.base
from typing import Union


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
    

    data_df = read_kfit_data_xlsx(join(ketchup_options['directory_data'],
                                       ketchup_options['filename_data']))

    data_dict = create_data_dict(data_df)
    
    # TODO: add code to set default values for secondary options if not in ketchup_options
 
    ketchup_model = create_initial_model(m_model, mech_df, data_dict, seedvalue=ketchup_options['seedvalue'],distribution=ketchup_options['distribution'],mech_type=ketchup_options['mechanism_type'])

    try:
        ketchup_model.name = str(ketchup_options['model_name'])
    except NameError:
        ketchup_model.name = 'unknown'
    
    experiments = []
    for i,d in enumerate(data_dict.keys()):
        ketchup_model.key = d
        ketchup_model.add_component(f'experiment{i}', Block(rule=ktools.core.create_sMB))
        experiments.append(eval(f"ketchup_model.experiment{i}"))
        print(f"dataset included: {d}")
    
    ketchup_model.add_component('block_list', pyomo.environ.Set(initialize = experiments))
    ketchup_model.reaction_rate = Constraint(ketchup_model.block_list,
                                                     ketchup_model.experiment0.REACTIONS,
                                                     rule=ktools.core.net_reaction_rate)
    
    try:
        basis_expt = str(ketchup_options['basis_id'])
    except KeyError:
        basis_expt = 'WT' if 'WT' in data_dict else data_dict[0]


    #ketchup_model.experiment0.vf.pprint()
    if ketchup_options['mechanism_type'] == 'elemental':
        def ref_conc_constraint(m,s):
            return m.experiment0.c[s] == 1.0
            
                #fix wild-type fluxes
        for r in data_dict[basis_expt]:
            ketchup_model.experiment0.rate[r].fix(data_dict[basis_expt][r][0])
             
        ketchup_model.ref_conc_constraint = Constraint(ketchup_model.experiment0.SPECIES, rule=ref_conc_constraint)        
        ketchup_model.obj = Objective( sense=pyomo.environ.minimize,
                     expr = sum(  b.error for b in experiments )               
                         )
        ketchup_model.vf_rate = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_F, rule = ktools.core.elemental_vf)
        ketchup_model.vr_rate = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_R, rule = ktools.core.elemental_vr)
        ketchup_model.es_net  = Constraint(ketchup_model.block_list, ketchup_model.ELEMENTALSTEP_F, rule = ktools.core.es_net_balance) 
     
    return ketchup_model


def solve_ketchup_model(ketchup_model: pyomo.core.base.PyomoModel.ConcreteModel,
                        ketchup_options: dict):
    """Returns solved pyomo Concrete Model using solver options file
    """
    
    solver = pyomo.environ.SolverFactory('ipopt')
    solver.options.option_file_name = ketchup_options['solver_opt_fn']

    try:
        results = solver.solve(ketchup_model, tee=True,symbolic_solver_labels=True);
        #total_error = sum( b.error.value for exps in experiments for b in exps[:] )
        status = str(results.Solver[0]['Termination condition']) 
    except ValueError:
        status = 'Status Error'

    return results




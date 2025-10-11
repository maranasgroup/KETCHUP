"""
FORmulate Kinetic model

Formulates the kinetic models for the solver
"""
from typing import Any

from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import pandas as pd
import re
import cobra.core.model

def create_data_dict(experiments_df: pd.DataFrame, data_selector: dict = None) -> dict:
    """
    Create data dictionary from the experiments in a form used in Pyomo.DAE computations
    For use with data for with stead state solutions.

    Parameters
    ----------
    experiments_df : pd.DataFrame
        DataFrame containing experimental data.
    data_selector: dict
        Dictionary by experiment ID that selects which subset of the data to use in the analysis.
        Currently, this feature is not yet implemented.

    Returns
    -------
    dict
        Nested dictionary of experimental data arranged by blocks.
    """
    experiment_keys = experiments_df['experiment ID'].unique().tolist()
    data = {}
    rev_data = {}

    for item in experiment_keys:#range(len(experiment_keys)):

        # set up the flux data (flx) as a dictionary
        _tmp = experiments_df[(experiments_df['experiment ID'] == item) & (experiments_df['type'] == 'flx')]

        t = {}

        for idx, row in _tmp.iterrows():
            t[row['rxn ID']] = (row['flx'], row['SD'])
        data[item] = t

        # set up rev items as second dictionary
        _tmp = experiments_df[(experiments_df['experiment ID'] == item) & (experiments_df['type'] == 'rev')]

        if not _tmp.empty:
            t = {}

            for idx, row in _tmp.iterrows():
                t[row['rxn ID']] = (row['flx'], row['SD'])
            rev_data[item+'_rev'] = t

    # For now, only return the 'flx' components. The 'rev' are not yet used.
    return data


def create_dynamic_data_dict(experiments_df: pd.DataFrame, data_selector: dict = None):
    """
    Create data dictionary from the experiments in a form used in Pyomo.DAE computations.
    For use with dynamic data that depend on time.

    Parameters
    ----------
    experiments_df : pd.DataFrame
        DataFrame containing experimental data.
    data_selector: dict
        Dictionary by experiment ID that selects which subset of the data to use in the analysis.
        Currently, this feature is not yet implemented.

    Returns
    -------
    dict
        Nested dictionary of experimental data arranged by blocks.
    """

    # As with static_data_dict, data blocks are indicated by experiment id

    # remove first two entries, since they are part of the meta-data
    from ktools.util import extract_from_brackets

    experiment_keys = experiments_df.loc[2:, 'experiment ID'].unique().tolist()

    # print (experiment_keys)

    data_dict = {}

    col_t0 = [col for col in experiments_df.columns if col.endswith("]_0")]

    col_e = [col for col in experiments_df.columns if experiments_df.loc[0, col] in (["e"])]  # dfr.iloc[0].isin(["e"])
    col_c = [col for col in experiments_df.columns if experiments_df.loc[0, col] in (["c"])]

    col_dependent = [col for col in experiments_df.columns if experiments_df.loc[1, col] in ["dependent"]]
    col_independent = [col for col in experiments_df.columns if experiments_df.loc[1, col] in ["independent"]]
    col_ignore = [col for col in experiments_df.columns if
                  experiments_df.loc[1, col] not in ["time", "independent", "dependent"]]
    col_time = [col for col in experiments_df.columns if experiments_df.loc[1, col] in ["time", "dependent"]]

    # TODO: reindex columns to ensure that Time will always come first

    # print (f"{col_dependent = }\n{col_ignore = }\n{col_e =}" )

    for k in experiment_keys:
        # set up the data as a dictionary of just the one experiment
        _df_tmp = experiments_df[(experiments_df['experiment ID'] == k)].reset_index(drop=True)
        # note: each slice now starts at row 0

        # filter all columns that are flagged to be ignored
        # TODO: add method to enable data selection at the individual experiment level.
        #       for now, the data_selector dictionary is not used
        _df_tmp = _df_tmp.drop(columns=col_ignore)

        # print (_df_tmp)
        _dict_t0_tmp = {}
        for item in col_t0:
            # print (extract_from_brackets(item))
            try:
                # print(extract_from_brackets(item), _tmp.loc[2,item])
                _dict_t0_tmp[extract_from_brackets(item)] = _df_tmp.loc[0, item]
            except Exception as e:
                #print(f"{e}")
                pass
            # print (_dict_t0_tmp)

        _dict_dep_tmp = {}
        # exract out the single time data column
        _dict_dep_tmp['Time'] = _df_tmp.loc[:, col_time[0]].values.tolist()
        #print (f"{col_dependent = }")
        for item in col_dependent:
            try:
                _dict_dep_tmp[extract_from_brackets(item)] = _df_tmp.loc[:, item].values.tolist()
            except Exception as e:
                print (f"{e}")
                pass
            #print (_dict_dep_tmp)

        data_dict.update({k: {"t0": _dict_t0_tmp, "time": _dict_dep_tmp}})

        # print (_tmp[col_dependent])

    return data_dict


def parse_time_delay(tde: int | float | list | dict, keys: list) -> dict:
    """Parses and validates the time delay input and converts it into a standardized dictionary.

    Parameters
    ----------
    tde:
        Time delay input. Can be a single value (applied to all experiments), or a
        list or dictionary (one per item).
    keys : list
        List of experiment IDs.

    Returns
    -------
    dict
        Dictionary of the time delays in a standard format using the experiment IDs as keys.

    Raises
    ------
    ValueError
        If 'tde' is not of a supported data type or is of the wrong length
    KeyError
        If provided experiment ID are missing from 'tde'.
    """

    dict_tde: dict[int, float] = {}
    if isinstance(tde, (int, float)):
        # make all the same
        for k in keys:
            dict_tde[k] = tde
    elif isinstance(tde, list):
        if len(tde) != len(keys):
            raise ValueError(f"'time_delay' option length does not agree with the number of experiments")
        else:
            for i,k in enumerate(keys):
                dict_tde[k] = tde[i]
    elif isinstance(tde, dict):
        for k in keys:
            if k not in tde:
                raise KeyError(f"Experiment ID '{k}' is not present in 'time_delay' option")
        for k in keys:
            dict_tde[k] = tde[k]
    else:
        raise ValueError(f"'time_delay' option is not of a correct data type: currently {type(tde)}")
    return dict_tde


### CONSTRAINT DEFINITIONS ###


#enz-sum
def enz_sum(m,r):
    """ Enzyme sums """
    rhs = 0
    for e in m.res[r][0]:
        rhs += m.e[e]
    return rhs == m.res[r][1]


#steady-state definition
def stoichiometry(m,s):
    """ steady-state definition """
    rhs = 0
    met = m.m_model.metabolites.get_by_id(s)
    for r in met.reactions:
        coeff = r.metabolites[met]
        rhs += m.rate[r.id]*coeff
    return rhs == 0.0


#elemental forward rate 
def elemental_vf(m,b,es):
    """ elemental forward rate """
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = m.kf[es]*b.e[df['reactant'].tolist()[0][0]]
    if len(df['reactant'].tolist()[0]) > 1:  rhs *= b.c[df['reactant'].tolist()[0][1]]
    return b.vf[es] == rhs


#elemental reverse rate
def elemental_vr(m,b,es):
    """ elemental reverse rate """
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = m.kr[es]*b.e[df['product'].tolist()[0][0]]
    if len(df['product'].tolist()[0]) > 1:  rhs *= b.c[df['product'].tolist()[0][1]]
    return b.vr[es] == rhs


#net reaction in elemental or MM
def net_reaction_rate(m,b,r):
    """ net reaction in elemental or MM """
    if m.mech_type == 'elemental':
        return Constraint.Skip
    elif m.mech_type == 'michaelis-menten':
        substrates = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'reactant'].tolist()
        products   = b.mech_df['product'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'product'].tolist()
        inh_ci  = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'competitive'].tolist()
        inh_uci = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'uncompetitive'].tolist()
        inh_nci = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'noncompetitive'].tolist()
        
        num = m.Kcat_f[r] * eval('*'.join([f"b.c[\'{c[0]}\']" for c in substrates])) / (eval('*'.join([f"m.KM_reactants[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).reactants]])))
        den = eval('+'.join( f"{abs(b.m_model.reactions.get_by_id(r).metabolites[b.m_model.metabolites.get_by_id(c)])}*b.c[\'{c}\']/m.KM_reactants[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).reactants]  ))
        #uncompetitive inhibition
        if inh_uci != [] or inh_nci != []:
            den = (den) * (1 + eval('+'.join(f"b.c[\'{c[0]}\']/m.KM_inhibitors[\'{r}+{c[0]}_uci\']" for c in inh_uci+inh_nci))  )  

        #competitive inhibition
        if inh_ci != [] or inh_nci != []:
            den += eval( '+'.join( [f"b.c[\'{c[0]}\']/m.KM_inhibitors[\'{r}+{c[0]}_ci\']" for c in inh_ci+inh_nci] ) )
                   
        if products != []:
            num -= m.Kcat_r[r] * eval('*'.join([f"b.c[\'{c[0]}\']" for c in products])) / (eval('*'.join([f"m.KM_products[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]])))
            den += eval('+'.join( f"{abs(b.m_model.reactions.get_by_id(r).metabolites[b.m_model.metabolites.get_by_id(c)])}*b.c[\'{c}\']/m.KM_products[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]  ))
        else:
            if r in [rxn.id for rxn in b.m_model.boundary if rxn.lower_bound < 0]: num *= -1    
        rhs = num * b.e[f"{r}_ENZ"]/den  
        
    return b.rate[r] == rhs


#elemental forward and reverse balance
def es_net_balance(m,b,es):
    """ elemental forward and reverse balance """
    check = es.split("_")[-1]
    if 'i' in check: return 0 == b.vf[es] - b.vr[es]
    rid = es[:-(len(check)+1)]
    return b.rate[rid] == b.vf[es] - b.vr[es]


### END CONSTRAINT DEFINITIONS ###

def create_initial_model(m_model: cobra.core.model.Model, mech_df : pd.DataFrame,
                         data: dict, seedvalue: int = None, time: float = None,
                         tde=0, k_thres: float = 10**5,
                         distribution: str = 'uniform', mech_type: str = 'elemental') -> pyomo.core.ConcreteModel:
    """
    Creates an initial Pyomo model for the problem to be solved.

    Parameters
    ----------
    m_model : cobra.core.model.Model
        The COBRApy model corresponding to the metabolic network.
    mech_df : pandas.DataFrame
        DataFrame containing mechanism information.
    data : dict
        Standardized data dictionary to use in the analysis.
    seedvalue : int, optional
        Seed for random number generation.
    time : list or float, optional
        List or float representing maximum value used for continuous time variable. Defaults to None.
    tde : float, optional
        Time delay. Defaults to 0.
    k_thres : float, optional
        Upper bound threshold for kinetic rate constants. Defaults to 10**5.
    distribution : str, optional
        Distribution type for initial random rate constants. Defaults to 'uniform'.
    mech_type : str, optional
        Mechanism type to follow for rate laws. Defaults to 'elemental'.

    Returns
    -------
    Pyomo.core.ConcreteModel
        The initialized Pyomo model.
    """

    if seedvalue == None:
        import warnings
        warnings.warn("Warning: seed value not provided. Setting one based on current time.")
        import time as timem
        seedvalue = int(timem.time() * 1000)  # ms since epoch
        # check validity of seed and assign value
        seedvalue = seedvalue % (2 ** 32 - 1)  # seed must be between 0 and 2**32 - 1
        print(f"Seed of {seedvalue} used.")
    np.random.seed(seedvalue)
    #generate random K values
    def randomize_Ks( mylist, mymodelvar, e_scale, distribution):
        plot_distribution = False
        dist_list = []
        for m in mylist:
            val = np.random.rand()
            if distribution == 'uniform' : val *= e_scale
            if distribution == 'log'     : val = 10**(val*np.log10(e_scale))
            mymodelvar[m] = val
            dist_list.append(val)
        if plot_distribution:
            import matplotlib.pyplot as plt
            plt.hist(dist_list,density=False,bins=30)
            plt.savefig('hist.png')
    model = ConcreteModel()
    model.mech_df = mech_df
    model.m_model = m_model
    model.data = data
    model.mech_type = mech_type

    e_scale = 5000
    ##INITIAL ELEMENTAL MODEL
    if mech_type == 'elemental':
        model.rxn_enz_sum = {rxn:[[item for sublist in mech_df['product'].loc[mech_df['rxn ID'] == rxn].tolist() for item in sublist if 'ENZ' in item],1] 
                   for rxn in [r.id for r in m_model.reactions]}
        #elemental steps
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        model.kf = Var( model.ELEMENTALSTEP_F, bounds=(0,k_thres) )
        model.kr = Var( model.ELEMENTALSTEP_R, bounds=(0,k_thres) )

        #distribution
        #e_scale = k_thres
        #randomize kf and kr
        
        randomize_Ks(elemental_steps, model.kf, e_scale, distribution)
        randomize_Ks(elemental_steps, model.kr, e_scale, distribution)
        
    elif mech_type == 'michaelis-menten':
        model.rxn_enz_sum = {rxn:[[f"{rxn}_ENZ"],1] for rxn in [r.id for r in m_model.reactions]}
        #MM formualation
        KMr = []; KMp = []; Ki = []
        inhib_dict = {"noncompetitive":"nci", "competitive":"ci", "uncompetitive":"uci"}
        for i,r in enumerate(mech_df['rxn ID'].tolist()):
            if mech_df.iloc[i]['type'] == 'reactant': KMr.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}')
            if mech_df.iloc[i]['type'] == 'product': KMp.append(f'{r}+{mech_df.iloc[i]["product"][0]}')
            if mech_df.iloc[i]['type'] in inhib_dict.keys(): Ki.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}_{inhib_dict.get(mech_df.iloc[i]["type"])}')
        KMi = [f"{i.split('_nci')[0]}_ci" for i in Ki if i.split('_')[-1] == 'nci']\
              +  [f"{i.split('_nci')[0]}_uci" for i in Ki if i.split('_')[-1] == 'nci']\
              +  [i for i in Ki if i.split('_')[-1] != 'nci']

        #Kinetic Parameter sets (reactants, products, inhibitors)
        model.KPr = Set( initialize = KMr )
        model.KPp = Set( initialize = KMp )
        model.KPi = Set( initialize = KMi )
        model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions])
        #variables
        model.KM_reactants  = Var( model.KPr, bounds = (0,k_thres) )
        model.KM_products   = Var( model.KPp, bounds = (0,k_thres) )
        model.KM_inhibitors = Var( model.KPi, bounds = (0,k_thres) )    
        model.Kcat_f        = Var( model.REACTIONS, bounds = (0,k_thres) )
        model.Kcat_r        = Var( model.REACTIONS, bounds = (0,k_thres) )
        model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )
        model.c = Var(model.SPECIES, bounds=(0,10**3), initialize=1 )        
        e_scale = k_thres

        #randomize KM
        randomize_Ks( KMr, model.KM_reactants , e_scale, distribution )
        randomize_Ks( KMp, model.KM_products  , e_scale, distribution )
        randomize_Ks( KMi, model.KM_inhibitors, e_scale, distribution )
    elif mech_type == 'custom':
        #iterate through mechanism file to find all rate laws and unique constants
        unique_consts = {'KCAT':[],'KM':[],'KI':[], 'KCONS':[],'rxns':{}}; tmp = []
        for i in mech_df.index:
            mech_info = mech_df.loc[i]['type']; rxn_name = mech_df.loc[i]['rxn ID']
            if mech_info not in ['reactant','product']:
                raw_rate_law = mech_info ; unique_consts['rxns'].update({rxn_name:raw_rate_law})
                for x in raw_rate_law.replace(' ', '').replace('*', '+').replace('/', '+').replace('(', '').replace(')','').replace('-', '+').split('+'):
                    if x not in tmp and not x.isdigit():
                #for x in raw_rate_law.replace(' ','').replace('*','+').replace('/','+').replace('(','').replace(')','').split('+'):
                #    if x not in tmp:
                        search_str = re.search(r'(?<=\[).+?(?=\])',x).group()
                        tmp.append(x); rxn_id = f"{rxn_name}+{search_str}"
                        if search_str == '[E]' : continue #ignore enzyme concentration
                        if 'KI' in x  : unique_consts['KI'].append(rxn_id)
                        if 'KCAT' in x: unique_consts['KCAT'].append(rxn_id)
                        if 'KM' in x  : unique_consts['KM'].append(rxn_id)    
                        if 'KCONS' in x  : unique_consts['KCONS'].append(rxn_id) 

        model.KCAT = Var(unique_consts['KCAT'],bounds=(0,k_thres))
        model.KM = Var(unique_consts['KM'],bounds=(0,k_thres))
        model.KI = Var(unique_consts['KI'],bounds=(0,k_thres))
        model.KCONS = Var(unique_consts['KCONS'],bounds=(0,k_thres))
        model.kinetic_parameters = {"kcat":model.KCAT , "Km":model.KM, "Ki":model.KI, "Kconsts":model.KCONS}
        randomize_Ks(unique_consts['KCAT'], model.KCAT, e_scale, distribution)
        randomize_Ks(unique_consts['KM'], model.KM, e_scale, distribution)
        randomize_Ks(unique_consts['KI'], model.KI, e_scale, distribution)
        randomize_Ks(unique_consts['KCONS'], model.KCONS, e_scale, distribution)
        model.unique_consts = unique_consts
        model.rxn_enz_sum = {rxn:[[f"{rxn}_ENZ"],1] for rxn in [r.id for r in m_model.reactions]}
        
    if time:
        try:
            maxTime = max(time)
            times = time
        except TypeError:
            maxTime = time
            times = [time]
    #time lag
    model.tde = tde
    #model.TIME = ContinuousSet( bounds=(0,maxTime), initialize=times )

    return model


## STATIC ##


#create static model block
def create_sMB(b):
    """ create static model block """
    data = b.model().data
    key = b.model().key
    mech_type = b.model().mech_type

    if key != 'WT':
        b.model().rxn_enz_sum[key][1] = 0

    model = create_sKM( b.model().rxn_enz_sum, b.model().m_model, b.model().mech_df, b.model().mech_type )
    
    if key != 'WT':
        b.model().rxn_enz_sum[key][1] = 1     
    
    model.error = Var(bounds=(0,None))  
    model.compute_error = Constraint(
         expr = model.error == sum(
            ( (model.rate[rxn]  - data[key][rxn][0] )**2/(data[key][rxn][1])**2 for rxn in data[key])
                                  )
    )
    
    return model    


def create_sKM(rxn_enz_sum, m_model, mech_df, mech_type, k_thres=10**10 ):
    """ create static kinetic model """
    model = ConcreteModel()
    model.mech_df = mech_df; model.res = rxn_enz_sum; model.m_model = m_model
    enz_names = [x for xs in list(rxn_enz_sum.values()) for x in xs[0]]
    model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )
    model.ENZYMES = Set (initialize = enz_names )    
    model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions] )
    model.c = Var(model.SPECIES, bounds=(0,10**3), initialize=1 )
    model.e = Var(model.ENZYMES, bounds=(0,1) )
    model.rate = Var(model.REACTIONS)
    
    if mech_type == 'elemental':
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        model.vf = Var( model.ELEMENTALSTEP_F, bounds=(0,None) )
        model.vr = Var( model.ELEMENTALSTEP_R, bounds=(0,None) )
    elif mech_type == 'michaelis-menten':
        pass        
        
    #CONSTRAINTS
    model.enz_sum = Constraint(model.REACTIONS, rule=enz_sum)  
    model.stoichiometry = Constraint(model.SPECIES,  rule=stoichiometry)        
    
    return model


## DYNAMIC ##


# create dynamic model block
def create_dMB(b):
    """ create dynamic model block """
    data = b.model().data
    key = b.model().key
    #print (f"{key = }")
    mech_type = b.model().mech_type
    tde = b.model().tde[key]
    #print(f"{tde = }")
    # set tde to 0 if found lag makes time less than 0
    #print (data[key]['time'])
    if data[key]['time']['Time'][0] - tde <  0:
        tde = data[key]['time']['Time'][0]
        print(f"#### experiment {key} results in time lag less than 0, resetting to tde to 0",flush=True)

    # set precision for rounding
    round_precision = 3
    #print ([t for t in data[key]['time']['Time'] ])
    time = [round(t-tde,round_precision) for t in data[key]['time']['Time'] ]
    # reconstruct data
    # TODO: allow more than one observation
    #print (list(data[key]['time'].keys()))
    error_key = list(data[key]['time'].keys())[1]
    #print (f"{error_key = }")
    error_data = [ (round(t-tde,round_precision),data[key]['time'][error_key][i]) \
                   for i,t in enumerate(data[key]['time']['Time']) ]
    #print (f"{data[key]= }")

    # insert (0,0) datapoint
    if tde != 0:
        time.insert(0,0)
        error_data.insert(0,(0,data[key]['t0'][error_key]))

    # setup enzyme sums
    # use gene_reaction_rules to associate data with the model. Enzyme ID values must match!
    for rxn in b.model().m_model.reactions:
        # right now only accepts singletons
        # TODO: update to allow multiple enzymes per reaction. Will require updates elsewhere
        pr = str(rxn.gene_reaction_rule)
        if pr:
            try:
                b.model().rxn_enz_sum[rxn.id][1] = data[key]['t0'][pr]
                #print (f"{rxn} catalyzed by {pr}, which at t0 = {data[key]['t0'][pr]}")
            except Exception as e:
                #print (e)
                print(f"{rxn} catalyzed by {pr} but no matching t0 value given. Check enzyme IDs. Assuming default enzyme concentration of 1")
                pass

    model = create_dKM(b.model().rxn_enz_sum, b.model().m_model, b.model().mech_df, b.model().mech_type, time,  data[key])
    model.data = data
    model.key = key
    model.tde = tde

    #reconfigure for time course
    model.error = Var(bounds=(0,None))
    # TODO: set up loop over multiple time data observations within single experiment
    model.compute_error = Constraint(
         expr = model.error == (1/len(error_data))*sum(
                       (model.c[t[0],error_key] - t[1])**2
                       for t in error_data
                                  )
    )
       
    model.del_component(model.TIME)
    return model


#create dynamic kinetic model
def create_dKM(rxn_enz_sum, m_model, mech_df, mech_type, time, data, k_thres=10**10 ):
    """ create dynamic kinetic model """
    model = ConcreteModel()
    model.mech_df = mech_df; model.res = rxn_enz_sum; model.m_model = m_model; model.data = data
    enz_names = [x for xs in list(rxn_enz_sum.values()) for x in xs[0]]
    model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )   
    model.ENZYMES = Set (initialize = enz_names ) 
    model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions] )    

    #KETCHUP_DYNAMIC additional variables
    maxtime = max(time)
    model.time = Set(initialize = time)
    model.TIME = ContinuousSet(bounds=(0,maxtime), initialize=time)
    model.e = Var(model.TIME, model.ENZYMES, bounds=(0,1) )    
    model.rate = Var(model.TIME, model.REACTIONS, bounds=(0,10*4))
    model.c = Var(model.TIME, model.SPECIES, bounds=(0,10**4), initialize=1 ) 
    model.dcdt = DerivativeVar(model.c, wrt=model.TIME)   
    
    if mech_type == 'elemental':
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        #KETCHUP_DYNAMIC added time variable to vf and vr
        model.vf = Var(model.TIME, model.ELEMENTALSTEP_F, bounds=(0,None) )
        model.vr = Var(model.TIME, model.ELEMENTALSTEP_R, bounds=(0,None) )

        model.enz_sum = Constraint(model.TIME, model.REACTIONS, rule=d_enz_sum)  
    elif mech_type == 'michaelis-menten':
        pass

    #CONSTRAINTS
    model.stoichiometry = Constraint(model.TIME, model.SPECIES,  rule=d_stoichiometry)
    
    model.enz_sum = Constraint(model.TIME, model.REACTIONS, rule=d_enz_sum)  
    
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(model,wrt=model.TIME,nfe=2*len(time),scheme="BACKWARD")
    
    return model


# TODO: define/move automatic way to transcribe written equation into symbolic form,
#       which is currently embedded in ketchup dynamic

# CUSTOM RATE LAWS
def custom_rate(b,t):
    """ create custom rates """
    reactions = [r.id for r in b.model().m_model.reactions]
    mech_df = b.model().mech_df
    return 0


#dynamic - enz-sum for each time-point
def d_enz_sum(m,t,r):
    """ dynamic enzyme sum. Covers each time-point """
    rhs = 0
    for e in m.res[r][0]:
        rhs += m.e[t,e]
    return rhs == m.res[r][1]


#dynamic - stoichiometric balance 
def d_stoichiometry(m,t,s):
    """ dynamic stoichiometric balance """
    rhs = 0
    met = m.m_model.metabolites.get_by_id(s)
    for r in met.reactions:
        coeff = r.metabolites[met]
        rhs += m.rate[t,r.id]*coeff
    return  m.dcdt[t,s] ==  rhs


#dynamic - MM form 
def d_net_reaction_rate(b,t,r):
    """ dynamic - michaelis-menten form """
    if b.model().mech_type == 'michaelis-menten':
        substrates = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'reactant'].tolist()
        products   = b.mech_df['product'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'product'].tolist()
        inh_ci  = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'competitive'].tolist()
        inh_uci = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'uncompetitive'].tolist()
        inh_nci = b.mech_df['reactant'].loc[b.mech_df['rxn ID'] == r].loc[b.mech_df['type'] == 'noncompetitive'].tolist()
        
        num = b.model().Kcat_f[r] * eval('*'.join([f"b.c[t,\'{c[0]}\']" for c in substrates])) / (eval('*'.join([f"b.model().KM_reactants[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).reactants]])))
        den = eval('*'.join( f"( 1 +{abs(b.m_model.reactions.get_by_id(r).metabolites[b.m_model.metabolites.get_by_id(c)])}*b.c[t,\'{c}\']/b.model().KM_reactants[\'{r}+{c}\'])" for c in [x.id for x in b.m_model.reactions.get_by_id(r).reactants]  ) )

        #uncompetitive inhibition
        if inh_uci != [] or inh_nci != []:
            den = (den) * (1 + eval('+'.join(f"b.c[t,\'{c[0]}\']/b.model().KM_inhibitors[\'{r}+{c[0]}_uci\']" for c in inh_uci+inh_nci))  )  

        #competitive inhibition
        if inh_ci != [] or inh_nci != []:
            den += eval( '+'.join( [f"b.c[t,\'{c[0]}\']/b.model().KM_inhibitors[\'{r}+{c[0]}_ci\']" for c in inh_ci+inh_nci] ) )
                   
        if products != []:
            num -= b.model().Kcat_r[r] * eval('*'.join([f"b.c[t,\'{c[0]}\']" for c in products])) / (eval('*'.join([f"b.model().KM_products[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]])))
            den += eval('*'.join( f"( 1 +{abs(b.m_model.reactions.get_by_id(r).metabolites[b.m_model.metabolites.get_by_id(c)])}*b.c[t,\'{c}\']/b.model().KM_products[\'{r}+{c}\'])" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]  ))
        else:
            if r in [rxn.id for rxn in b.m_model.boundary if rxn.lower_bound < 0]: num *= -1    
        rhs = num * b.e[t,f"{r}_ENZ"]/(den - 1)
    else:
        return Constraint.Skip
        
    return b.rate[t,r] == rhs  

"""
FORmulate Kinetic model

Formulates the kinetic models for the solver
"""

from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import pandas as pd

def create_data_dict(experiments_df):
    """ Create data dictionary from the experiments in a form used in Pyomo.DAE computations
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
    return (data)
###CONSTRAINT DEFINITIONS ###
#enz-sum
def enz_sum(m,r):
    rhs = 0 
    for e in m.res[r][0]:
        rhs += m.e[e]
    return rhs == m.res[r][1]

#steady-state definition
def stoichiometry(m,s):
    rhs = 0
    met = m.m_model.metabolites.get_by_id(s)
    for r in met.reactions:
        coeff = r.metabolites[met]
        rhs += m.rate[r.id]*coeff
    return rhs == 0.0
#elemental forward rate 
def elemental_vf(m,b,es):
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = m.kf[es]*b.e[df['reactant'].tolist()[0][0]]
    if len(df['reactant'].tolist()[0]) > 1:  rhs *= b.c[df['reactant'].tolist()[0][1]]
    return b.vf[es] == rhs
#elemental reverse rate
def elemental_vr(m,b,es):
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = m.kr[es]*b.e[df['product'].tolist()[0][0]]
    if len(df['product'].tolist()[0]) > 1:  rhs *= b.c[df['product'].tolist()[0][1]]
    return b.vr[es] == rhs
#net reaction in elemental or MM
def net_reaction_rate(m,b,r):
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

def es_net_balance(m,b,es):
    check = es.split("_")[-1]
    if 'i' in check: return 0 == b.vf[es] - b.vr[es]
    rid = es[:-(len(check)+1)]
    return b.rate[rid] == (b.vf[es] - b.vr[es]) * m.rxn_enz_sum[rid][1]
    
def create_initial_model(m_model,mech_df,data,seedvalue=None,time=None,k_thres=10**10,mech_type='elemental'):
    if seedvalue == None:
        import warnings
        warnings.warn("Warning: seed value not provided. Setting one based on current time.")
        from datetime import datetime
        seedvalue = 0
        current_datetime = datetime.now()
        timestamp = int(round(current_datetime.timestamp()))
        seedvalue = timestamp
        print(f"Seed of {seedvalue} used.")
    np.random.seed(seedvalue)
    #generate random K values
    def randomize_Ks( mylist, mymodelvar, e_scale):
        for m in mylist:
            val = np.random.random() * e_scale
            mymodelvar[m] = val

    model = ConcreteModel()
    model.mech_df = mech_df
    model.m_model = m_model
    model.data = data
    model.mech_type = mech_type
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

        e_scale = k_thres
        #randomize kf and kr
        for es in elemental_steps:
            valf = np.random.random() * e_scale
            valr = np.random.random() * e_scale
            #print(f'{e_step} = {val}')
            model.kf[es] = valf
            model.kr[es] = valr   
        randomize_Ks(elemental_steps, model.kf, e_scale)
        randomize_Ks(elemental_steps, model.kr, e_scale)
        
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
        model.c = Var(model.SPECIES, bounds=(0,1000), initialize=1 )        
        e_scale = k_thres

        #randomize KM
        randomize_Ks( KMr, model.KM_reactants , e_scale )
        randomize_Ks( KMp, model.KM_products  , e_scale )
        randomize_Ks( KMi, model.KM_inhibitors, e_scale )
    if time:
        try:
            maxTime = max(time)
            times = time
        except TypeError:
            maxTime = time
            times = [time]
        model.TIME = ContinuousSet( bounds=(0,maxTime), initialize=times )
    return model

#create static model block

def create_sMB(b):
    data = b.model().data
    key = b.model().key
    mech_type = b.model().mech_type
    
    p_dict = {'LDH_L':{'LDH_L':0},
              'ALCD2x;ACALDRXN':{'ALCD2x':0.5,'ACALDRXN':0.5},
              'ALCD2x;ACALDRXN;LDH_L':{'ALCD2x':0.5,'ACALDRXN':0.5,'LDH_L':0},
              'POR2_i':{'POR2_i':0},
              'PTAr':{'PTAr':0},
              'LDH_L;PTAr':{'LDH_L':0,'PTAr':0},
              'NADOX':{'NADOX':0},
              'NADOX;LDH_L;PFL;PTAr;ACKr':{'LDH_L':0,'PFL':0,'ACKr':0,'NADOX':0.5,'PTAr':0},
              'EXCH_h2(e)':{'EXCH_h2(e)':0,'NADOX':0.1},                  
              'ALCD2x;ACALDRXN;EXCH_lac-l(e)':{'ALCD2x':0.2,'ACALDRXN':0.2,'LDH_L':5}
              }
              
    if 'WT' not in key:
        for d in p_dict[key]:
            b.model().rxn_enz_sum[d][1] = p_dict[key][d]

    model = create_sKM( b.model().rxn_enz_sum, b.model().m_model, b.model().mech_df, b.model().mech_type )
     
    model.error = Var(bounds=(0,None))
    
    model.rate['EXCH_cellobiose(e)'].bounds = (None,-0.1)
    #yield calculation - need to program objective function can be toggled for yield or rates 
    model.y = Var(model.REACTIONS)  
    def yield_calc(m,r):
        return m.y[r] == 100 * m.rate[r] / (-m.rate['EXCH_cellobiose(e)'])
    model.yields = Constraint(model.REACTIONS, rule=yield_calc)
    
    model.compute_error = Constraint(
         expr = model.error == sum(
            ( (model.y[rxn]  - data[key][rxn][0] )**2/(data[key][rxn][1])**2 for rxn in data[key])
                                  )
    )
    
    
    if 'WT' not in key:
        for d in p_dict[key]:
            b.model().rxn_enz_sum[d][1] = 1 
            
            
    return model    

def create_sKM( rxn_enz_sum, m_model, mech_df, mech_type, k_thres=10**10 ):
    model = ConcreteModel()
    model.mech_df = mech_df; model.res = rxn_enz_sum; model.m_model = m_model; model.rxn_enz_sum = rxn_enz_sum
    enz_names = [x for xs in list(rxn_enz_sum.values()) for x in xs[0]]
    model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )
    model.ENZYMES = Set (initialize = enz_names )    
    model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions] )
    model.c = Var(model.SPECIES, bounds=(0,1000), initialize=1 )
    model.e = Var(model.ENZYMES, bounds=(0,10) )
    model.rate = Var(model.REACTIONS)
    
    if mech_type == 'elemental':
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        model.vf = Var( model.ELEMENTALSTEP_F, bounds=(0,None) )
        model.vr = Var( model.ELEMENTALSTEP_R, bounds=(0,None) )
                    
    elif mech_type == 'michaelis-menten':
        KMr = []; KMp = []; Ki = []
        inhib_dict = {"noncompetitive":"nci", "competitive":"ci", "uncompetitive":"uci"}
        for i,r in enumerate(mech_df['rxn ID'].tolist()):
            if mech_df.iloc[i]['type'] == 'reactant': KMr.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}')
            if mech_df.iloc[i]['type'] == 'product': KMp.append(f'{r}+{mech_df.iloc[i]["product"][0]}')
            if mech_df.iloc[i]['type'] in inhib_dict.keys(): Ki.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}_{inhib_dict.get(mech_df.iloc[i]["type"])}')
        KMi = [f"{i.split('_nci')[0]}_ci" for i in Ki if i.split('_')[-1] == 'nci']\
          +  [f"{i.split('_nci')[0]}_uci" for i in Ki if i.split('_')[-1] == 'nci']\
          +  [i for i in Ki if i.split('_')[-1] != 'nci']
        
        model.KPr = Set( initialize = KMr )
        model.KPp = Set( initialize = KMp )
        model.KPi = Set( initialize = KMi )
        #model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions])
        #variables
        model.KM_reactants  = Var( model.KPr, bounds = (0,k_thres) )
        model.KM_products   = Var( model.KPp, bounds = (0,k_thres) )
        model.KM_inhibitors = Var( model.KPi, bounds = (0,k_thres) )    
        model.Kcat_f        = Var( model.REACTIONS, bounds = (0,k_thres) )
        model.Kcat_r        = Var( model.REACTIONS, bounds = (0,k_thres) )         
        

    #CONSTRAINTS
    model.enz_sum = Constraint(model.REACTIONS, rule=enz_sum)  
    model.stoichiometry = Constraint(model.SPECIES,  rule=stoichiometry)        
    
    return model




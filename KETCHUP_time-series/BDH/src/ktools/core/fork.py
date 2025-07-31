"""
FORmulate Kinetic model

Formulates the kinetic models for the solver
"""

from pyomo.environ import *
from pyomo.dae import *
import numpy as np
import pandas as pd
eps = 1e-8
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
    elif r != 'BDH':
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
            #num -= m.Kcat_r[r] * eval('*'.join([f"b.c[\'{c[0]}\']" for c in products])) / (eval('*'.join([f"m.KM_products[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]])))
            den += eval('+'.join( f"{abs(b.m_model.reactions.get_by_id(r).metabolites[b.m_model.metabolites.get_by_id(c)])}*b.c[\'{c}\']/m.KM_products[\'{r}+{c}\']" for c in [x.id for x in b.m_model.reactions.get_by_id(r).products]  ))
        else:
            if r in [rxn.id for rxn in b.m_model.boundary if rxn.lower_bound < 0]: num *= -1    
        rhs = num * b.e[f"{r}_ENZ"]/den  
        
    return b.rate[r] == rhs

def es_net_balance(m,b,es):
    check = es.split("_")[-1]
    if 'i' in check: return 0 == b.vf[es] - b.vr[es]
    rid = es[:-(len(check)+1)]
    return b.rate[rid] == b.vf[es] - b.vr[es]
    
def create_initial_model(m_model,mech_df,data,seedvalue=None,time=None,tde = 0, k_thres=10**5,distribution='uniform',mech_type='elemental'):
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

    e_scale = k_thres
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
        model.KM_reactants  = Var( model.KPr, bounds = (eps,k_thres) )
        model.KM_products   = Var( model.KPp, bounds = (eps,k_thres) )
        model.KM_inhibitors = Var( model.KPi, bounds = (eps,k_thres) )    
        model.Kcat_f        = Var( model.REACTIONS, bounds = (eps,k_thres) )
        model.Kcat_r        = Var( model.REACTIONS, bounds = (eps,k_thres) )
        model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )
        model.c = Var(model.SPECIES, bounds=(0,10**3), initialize=1 )        
        e_scale = k_thres

        #randomize KM
        randomize_Ks( KMr, model.KM_reactants , e_scale, distribution )
        randomize_Ks( KMp, model.KM_products  , e_scale, distribution )
        randomize_Ks( KMi, model.KM_inhibitors, e_scale, distribution )

        model.K_decomp = Var(bounds=(0,k_thres))

        
        model.kinetic_parameters = {"kcat_f":model.Kcat_f ,"kcat_r":model.Kcat_r,
                                     "Km_r" :model.KM_reactants , "Ki" :model.KM_inhibitors, "Km_p":model.KM_products,
                                     "K_decomp" : model.K_decomp
                                    }
    elif mech_type == 'custom':
        #define rate equations here
        #kinetic parameters should be split into 3 categories:
        #1 a maximum kinetic constant eg kcat
        #2 kinetic parameters capturing substrate/product interaction
        #3 kinetic parameters capturing regulatory interaction
        model.kcat = Var(bounds = (0,k_thres) )
        #define Km 
        model.reactants = Set(initialize = ["nad","hco2"])
        model.Km_r   =  Var(model.reactants, bounds = (0,k_thres) )
        #define Ki
        model.inhibitors = Set(initialize = ['nadh'])
        model.Ki     =  Var(model.inhibitors, bounds=(0,k_thres) )
        model.K_decomp = Var(bounds=(0,k_thres))
        
        #grouping the kinetic parameters in a list
        model.kinetic_parameters = {"kcat":model.kcat , "Km_r" : model.Km_r , "Ki" : model.Ki , "K_decomp" : model.K_decomp}
        #only for this occasion
        model.rxn_enz_sum = {rxn:[[f"{rxn}_ENZ"],1] for rxn in [r.id for r in m_model.reactions]}
        randomize_Ks(["nad","hco2"], model.Km_r, e_scale, distribution)
        randomize_Ks(['nadh'], model.Ki, e_scale, distribution)
        
        
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
    
    
#create static model block
def create_sMB(b):
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

def create_sKM( rxn_enz_sum, m_model, mech_df, mech_type, k_thres=10**10 ):
    model = ConcreteModel()
    model.mech_df = mech_df; model.res = rxn_enz_sum; model.m_model = m_model
    enz_names = [x for xs in list(rxn_enz_sum.values()) for x in xs[0]]
    model.SPECIES = Set( initialize = [m.id for m in m_model.metabolites] )
    model.ENZYMES = Set (initialize = enz_names )    
    model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions] )
    model.c = Var(model.SPECIES, bounds=(0,10**3), initialize=1 )
    model.e = Var(model.ENZYMES, bounds=(0,100) ) #enzyme total should be unbounded for upper-bound, we set 100 to be an unrealistic overexpression value facilitate solving
    model.rate = Var(model.REACTIONS)
    
    if mech_type == 'elemental':
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        model.vf = Var( model.ELEMENTALSTEP_F, bounds=(0,None) )
        model.vr = Var( model.ELEMENTALSTEP_R, bounds=(0,None) )
    elif mech_type == 'michaelis-menten':
        pass
        #NOTE 08/11/2024 this section should be in the create intial model function because kinetic parameters are consistent throughout all blocks
        #KMr = []; KMp = []; Ki = []
        #inhib_dict = {"noncompetitive":"nci", "competitive":"ci", "uncompetitive":"uci"}
        #for i,r in enumerate(mech_df['rxn ID'].tolist()):
        #    if mech_df.iloc[i]['type'] == 'reactant': KMr.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}')
        #    if mech_df.iloc[i]['type'] == 'product': KMp.append(f'{r}+{mech_df.iloc[i]["product"][0]}')
        #    if mech_df.iloc[i]['type'] in inhib_dict.keys(): Ki.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}_{inhib_dict.get(mech_df.iloc[i]["type"])}')
        #KMi = [f"{i.split('_nci')[0]}_ci" for i in Ki if i.split('_')[-1] == 'nci']\
        #  +  [f"{i.split('_nci')[0]}_uci" for i in Ki if i.split('_')[-1] == 'nci']\
        #  +  [i for i in Ki if i.split('_')[-1] != 'nci']
        
        #model.KPr = Set( initialize = KMr )
        #model.KPp = Set( initialize = KMp )
        #model.KPi = Set( initialize = KMi )
        #model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions])
        #variables
        #model.KM_reactants  = Var( model.KPr, bounds = (0,k_thres) )
        #model.KM_products   = Var( model.KPp, bounds = (0,k_thres) )
        #model.KM_inhibitors = Var( model.KPi, bounds = (0,k_thres) )    
        #model.Kcat_f        = Var( model.REACTIONS, bounds = (0,k_thres) )
        #model.Kcat_r        = Var( model.REACTIONS, bounds = (0,k_thres) )         
        
        
    #CONSTRAINTS
    model.enz_sum = Constraint(model.REACTIONS, rule=enz_sum)  
    model.stoichiometry = Constraint(model.SPECIES,  rule=stoichiometry)        
    
    return model


#KETCHUP_DYNAMIC create dynamic model block 
def create_dMB(b):
    data = b.model().data
    key = b.model().key
    mech_type = b.model().mech_type
    #lag time 09/06/2024
    tde = b.model().tde
    print(f'tde = {tde}')
    #set tde to 0 if found lag makes time less than 0
    if data[key]['time'][0][0] - tde <  0: 
        tde = data[key]['time'][0][0]
        print(f"#### experiment {key} has time lag less than 0 , resetting to tde to 0",flush=True) 
    
    time = [round(t[0]-tde,3) for t in data[key]['time'] ]
    #reconstruct data
    error_data = [ (round(t[0]-tde,3),t[1]) for t in data[key]['time'] ]
    
    #insert (0,0) datapoint
    if tde != 0:
        time.insert(0,0)
        error_data.insert(0,(0,data[key]['nadh']))

    
    
    ##FOR BDH
    b.model().rxn_enz_sum['BDH'][1] = 1 if data[key]['enzyme'] > 1 else 0.1  

    model = create_dKM( b.model().rxn_enz_sum, b.model().m_model, b.model().mech_df, b.model().mech_type, time,  data[key])
    model.data = data
    model.key = key
    model.tde = tde
    model.c[0,'nad'] = data[key]['nad']
    
    model.c[0,'nadh'] = data[key]['nadh']
    model.c[0,'nadh'].fixed = True    
    
    model.c[0,'acetoin'] = data[key]['acetoin']
    model.c[0,'acetoin'].fixed = True

    model.c[0,'23bd'] = data[key]['23bd']
    model.c[0,'23bd'].fixed = True
    
    #reconfigure for time course
    model.error = Var(bounds=(0,None))  
    model.compute_error = Constraint(
         expr = model.error == sum(
                       (model.c[t[0],'nadh'] - t[1])**2
                       for t in error_data
                                  )
    )
       
    model.del_component(model.TIME)
    return model
    
#KETCHUP_DYNAMIC create dynamic kinetic model
def create_dKM( rxn_enz_sum, m_model, mech_df, mech_type, time, data, k_thres=10**10 ):

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
    model.dcdt = DerivativeVar(model.c, wrt=model.TIME, within=Reals)   
    
    if mech_type == 'elemental':
        elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
        model.ELEMENTALSTEP_F = Set( initialize = elemental_steps )
        model.ELEMENTALSTEP_R = Set( initialize = elemental_steps )
        #KETCHUP_DYNAMIC added time variable to vf and vr
        model.vf = Var(model.TIME, model.ELEMENTALSTEP_F, bounds=(0,None) )
        model.vr = Var(model.TIME, model.ELEMENTALSTEP_R, bounds=(0,None) )

        
        model.enz_sum = Constraint(model.TIME, model.REACTIONS, rule=d_enz_sum)  
    elif mech_type == 'michaelis-menten':
        model.enz_sum = Constraint(model.TIME, model.REACTIONS, rule=d_enz_sum)  
        pass
        #NOTE 08/11/2024 this section should be in the create intial model function because kinetic parameters are consistent throughout all blocks
        #KMr = []; KMp = []; Ki = []
        #inhib_dict = {"noncompetitive":"nci", "competitive":"ci", "uncompetitive":"uci"}
        #for i,r in enumerate(mech_df['rxn ID'].tolist()):
        #    if mech_df.iloc[i]['type'] == 'reactant': KMr.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}')
        #    if mech_df.iloc[i]['type'] == 'product': KMp.append(f'{r}+{mech_df.iloc[i]["product"][0]}')
        #    if mech_df.iloc[i]['type'] in inhib_dict.keys(): Ki.append(f'{r}+{mech_df.iloc[i]["reactant"][0]}_{inhib_dict.get(mech_df.iloc[i]["type"])}')
        #KMi = [f"{i.split('_nci')[0]}_ci" for i in Ki if i.split('_')[-1] == 'nci']\
        #  +  [f"{i.split('_nci')[0]}_uci" for i in Ki if i.split('_')[-1] == 'nci']\
        #  +  [i for i in Ki if i.split('_')[-1] != 'nci']
       # 
       # model.KPr = Set( initialize = KMr )
       # model.KPp = Set( initialize = KMp )
       # model.KPi = Set( initialize = KMi )
       # #model.REACTIONS = Set( initialize = [r.id for r in m_model.reactions])
       # #variables
       # model.KM_reactants  = Var( model.KPr, bounds = (0,k_thres) )
       # model.KM_products   = Var( model.KPp, bounds = (0,k_thres) )
       # model.KM_inhibitors = Var( model.KPi, bounds = (0,k_thres) )    
       # model.Kcat_f        = Var( model.REACTIONS, bounds = (0,k_thres) )
       # model.Kcat_r        = Var( model.REACTIONS, bounds = (0,k_thres) )    
             
    #08/05/2024 CUSTOM RATE LAWS
    #NOTE: We will have to build a function to aggregate related kinetic parameters depending on type
    #might have to build a dictionary to "link" each custom kinetic parmeters

    #CONSTRAINTS
    model.stoichiometry = Constraint(model.TIME, model.SPECIES,  rule=d_stoichiometry)
    
    model.data = data
    def nad_conc_bounds(m,t):
        return (m.data['nad']*1 + m.data['nadh']) - m.c[t,'nad'] >= 0 
    def nadh_conc_bounds(m,t):
        return (m.data['nad']*1 + m.data['nadh']) - m.c[t,'nadh'] >= 0 
        
    def acetoin_conc_bounds(m,t):
        return (m.data['acetoin']*1 + m.data['23bd']) - m.c[t,'acetoin'] >= 0
    def bd23_conc_bounds(m,t):
        return (m.data['acetoin']*1 + model.data['23bd']) - m.c[t,'23bd'] >= 0

               
    model.nad_conc = Constraint(model.TIME, rule=nad_conc_bounds)
    model.nadh_conc = Constraint(model.TIME, rule=nadh_conc_bounds)
    model.acetoin_conc = Constraint(model.TIME, rule=acetoin_conc_bounds)    
    model.bd23_conc = Constraint(model.TIME, rule=bd23_conc_bounds)  
    
    #for dynamic data
    #colloc = TransformationFactory('dae.collocation')  
    #colloc.apply_to(model,ncp = 3,nfe = (len(time)), wrt=model.TIME, scheme='LAGRANGE-RADAU') 
    
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(model,wrt=model.TIME,nfe=2*(len(time)),scheme='BACKWARD')
    
    return model
def fdh_rate(b,t):
    #enzyme_conc = 22 / 1000   #units are in mM, as stated in Schmidt enzyme is 22 uM diluted in 1:20
    #enzyme_conc = 1
    #enzyme_conc = 2.514 * 10 ** -4
    num = b.model().kcat * b.c[t,'nad'] * b.c[t,'hco2'] * b.e[t,'FDH_ENZ'] 
    den = b.model().Km_r['nad'] * b.model().Km_r['hco2'] + \
          b.model().Km_r['hco2'] * b.c[t,'nad'] + \
          b.model().Km_r['nad'] * b.c[t,'hco2'] + \
          b.c[t,'nad'] * b.c[t,'hco2'] + \
          (b.model().Km_r['nad']*b.model().Km_r['hco2'])/(b.model().Ki['nadh'])*b.c[t,'nadh'] + \
          (b.model().Km_r['nad']/b.model().Ki['nadh']) * b.c[t,'hco2'] * b.c[t,'nadh']
    return b.rate[t,'FDH'] == num/den    
def nadh_decomp(b,t):
    return b.rate[t,'NADH_decomp'] == b.model().K_decomp * b.c[t,'nadh']
    
#enz-sum
def d_enz_sum(m,t,r):
    rhs = 0 
    for e in m.res[r][0]:
        rhs += m.e[t,e]
    return rhs == m.res[r][1]

#steady-state definition
def d_stoichiometry(m,t,s):
    rhs = 0
    met = m.m_model.metabolites.get_by_id(s)
    for r in met.reactions:
        coeff = r.metabolites[met]
        rhs += m.rate[t,r.id]*coeff
    return  m.dcdt[t,s] ==  rhs
    
#elemental forward rate 
def d_elemental_vf(b,t,es):
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = b.model().kf[es]*b.e[t,df['reactant'].tolist()[0][0]]
    if len(df['reactant'].tolist()[0]) > 1:  rhs *= b.c[t,df['reactant'].tolist()[0][1]]
    return b.vf[t,es] == rhs
#elemental reverse rate
def d_elemental_vr(b,t,es):
    s = es.split('_')[-1]
    step = es[:-(len(s)+1)]
    df = b.mech_df.loc[(b.mech_df['rxn ID']==step)&(b.mech_df['step ID']==s)]
    rhs = b.model().kr[es]*b.e[t,df['product'].tolist()[0][0]]
    if len(df['product'].tolist()[0]) > 1:  rhs *= b.c[t,df['product'].tolist()[0][1]]
    return b.vr[t,es] == rhs  
      
def d_es_net_balance(b,t,es):
    check = es.split("_")[-1]
    if 'i' in check: return 0 == b.vf[t,es] - b.vr[t,es]
    rid = es[:-(len(check)+1)]
    return b.rate[t,rid] == b.vf[t,es] - b.vr[t,es]   

#net reaction in elemental or MM
def d_net_reaction_rate(b,t,r):

    if r != 'BDH':
        return Constraint.Skip
    elif b.model().mech_type == 'michaelis-menten':
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
        
    return b.rate[t,r] == rhs     

"""
Parsing of dataframes for import.

Parsing of K-FIT derived dataframes of model and mechanism data.

"""
import cobra
import pandas as pd

def parse_kfit_model_df(df_metabolite, df_reaction, model_name='model', debug=False):
    """ Parse the DataFrames of the K-FIT model inputs and convert into COBRApy model. 
    """

    m = cobra.Model(model_name)

    # serialize the metabolites, using (temporary) names that we know will be acceptable python variable names
    tmp_counter = 1
    for idx, row in df_metabolite.iterrows():
        item = row['ID']
        if item not in m.metabolites:
            met_var = f"M_{tmp_counter}"
            # note: K-FIT spreadsheets do not track charges
            vars()[met_var] = cobra.Metabolite(item,formula=row['Formula'],name=row['Name'])
            if row['Name'] is None:
                vars()[met_var].name = row['ID']
            # compartment hack for simple systems
            # TODO: make general for more complex systems
            if '_e' in item[-2:] or '[e]' in item [-3:]:
                comp='e'
            elif '_m' in item[-2:]:
                comp='m'
            else:
                comp='c'
            vars()[met_var].compartment=comp
            # TODO: set sbo for exchange metabolites
            vars()[met_var].annotation={'sbo':'SBO:0000247'}
            m.add_metabolites([vars()[met_var]])
            tmp_counter += 1

    #display(model)

    # serialize the reactions, using (temporary) names that we know will be acceptable python variable names
    tmp_counter = 1
    for idx, row in df_reaction.iterrows():
        item = str(row['Rxn ID'])
        if item not in m.reactions:
            rxn_var = f"R_{tmp_counter}"
            vars()[rxn_var] = cobra.Reaction(item, name=row['Rxn name'])
            m.add_reactions([vars()[rxn_var]])
            vars()[rxn_var].reaction = (str(row['Rxn Formula']))
            #set upper and lower bounds. done after parsing the formula as the arrows there will affect 
            #  the bounds and K-FIT input files are not consistent on arrow types
            vars()[rxn_var].lower_bound=row['Lower bound']
            vars()[rxn_var].upper_bound=row['Upper bound']
            # TODO: add sbo to reactions
            tmp_counter += 1

    if debug:
        print(f'COBRA model with name {model_name} created')
    return m

def parse_kfit_mech_df(df_mechanism, m_model=None, model_name='model',
                       mech_type='elemental', debug=False, testing_flag=False):
    """ Parse the DataFrames of the K-FIT mechanism inputs and convert into COBRApy model.
    Uses the COBRA model so validate the reaction and metabolite ids.
    """
    import cobra

    # first validate ID
    if m_model:
        rxn_ids = df_mechanism['ID'].squeeze().tolist()
        m_model_rxn_ids = [x.id for x in m_model.reactions]
        erronious_rxn_type1 = [x for x in rxn_ids if x not in m_model_rxn_ids]
        if erronious_rxn_type1:
            print (f"Mechanism has reaction id(s) inconsistent with the model file: {erronious_rxn_type1}")
        erronious_rxn_type2 = [x for x in m_model_rxn_ids if x not in rxn_ids]
        if erronious_rxn_type2:
            print (f"Mechanism file is missing reaction id from the model file: {erronious_rxn_type2}")

    # process in the same order as the COBRA reactions, if present
    if m_model:
        rxn_id_list = [x.id for x in m_model.reactions]
    else:
        rxn_id_list = df_mechanism['ID'].squeeze().tolist()
        
    def df_join_reactions(df=pd.DataFrame(), data=None):
        """ Add items to the dataframe or initialize the data frame if empty
        
        """
        chart_index = ['rxn ID','step ID','reactant','product','mechanism','type']
        if df.empty and data == None:
            df = pd.DataFrame(columns = chart_index)
        else:
            df = pd.DataFrame(data, columns = chart_index)
        return df

    def enzyme_complex_name(name, cmp):
        """ Convert list based bindings into readable enzyme complexes. Enables 
            having informative names instead of generic entries like R1_enz_complex_1
        """

        output = str(name) # TODO: update using join
        if cmp:
            for item in cmp:
                output = output + "+" + str(item)
        return output
    
    def generate_elemental_reactions(rxn_df, base_enzyme_name=None):
        mechanism = rxn_df.mechanism
        base_rxn_name = rxn_df.ID
        
        if not base_enzyme_name:
            # note that this step assumes that enzymes are independent for each reaction.
            # if some are the same, constraints will need to be added elsewhere
            base_enzyme_name = f"{rxn_df.ID}_ENZ"

        # TODO: validate all metabolites included

        rxn_list = []

        df_new_reactions = df_join_reactions()
        
        # substrate binding order
        sbo = rxn_df['SBO']
        if sbo:
            sbo = sbo.split(';')
        # product release order
        pro = rxn_df['PRO']
        if pro:
            pro = pro.split(';')
        if debug:
            print (sbo,pro)

        enz = base_enzyme_name

        # reaction steps
        
        if mechanism.lower() == 'seq':
            # steps are "substrate binding(s), catalytic, product release(s)"
            #
            # e.g., A + B -> P + Q via E with sbo: A, B & pro: P, Q
            #   _0 E + A -> EA : binding
            #   _1 EA + B -> EAB : binding
            #   _2 EAB -> EPQ : catalytic
            #   _3 EPQ -> EQ + P : release
            #   _4 EQ -> E + Q : release

            # note that constraints on exchange reactions are NOT set here. the full set of elemental
            #    reactions are returned and then set to 0, 1, etc. when constructing the kinetic model
            #    parameterization
            
            # substrate binding
            if sbo:
                # note: the order adds to the right
                for count,item in enumerate(sbo):
                    current_react_list = [enz,item]
                    enz = enzyme_complex_name(base_enzyme_name, sbo[:count+1])
                    current_prod_list  = [enz]
                    current_rxn_name = str(count)
                    
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'binding'])
                    #print (current_rxn_name, current_react_list, current_prod_list)
            # catalytic
            count = len(sbo)
            current_react_list = [enz]
            enz = enzyme_complex_name(base_enzyme_name, pro)
            current_prod_list = [enz]
            current_rxn_name = str(count)
            
            # add new data frame row
            rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'catalytic'])
            #print (current_rxn_name, current_react_list, current_prod_list)

            # product release
            if pro:
                # note: the order on the right pops off from left, which is the
                #   opposite of the add to the right for binding
                count_start = len(sbo)+1
                for count,item in enumerate(pro):
                    current_react_list = [enz]
                    enz = enzyme_complex_name(base_enzyme_name, pro[count+1:])
                    current_prod_list = [enz,item]
                    current_rxn_name = str(count+count_start)
                    
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'release'])
                    #print (current_rxn_name, current_react_list, current_prod_list)
        
        elif mechanism.lower() == 'ppg':
        #if mechanism.lower() == 'seq': # for quick testing
            # ping pong alternates from reactant to product and the len of sbo and pro must be the same

            # e.g. A + B -> P + Q via E with sbo: A, B & pro: P, Q
            # _0 E + A -> EA : binding
            # _1 EA -> EP : catalytic
            # _2 EP -> E1 + P : release
            # _3 E1 + B -> E1B : binding
            # _4 E1B -> E1Q : catalytic
            # _5 E1Q -> E + Q : release

            if len(sbo) != len(pro):
                print (f"Warning in ping-pong reaction {base_rxn_name}: substrates and product counts not equal.")
            else:

                for count,item in enumerate(sbo):
                    if count > 0:
                        current_base_enzyme_name = f"{base_enzyme_name}_{count}"
                    else:
                        current_base_enzyme_name = base_enzyme_name

                    # substrate binding
                    current_react_list = [enz,sbo[count]]
                    enz = enzyme_complex_name(current_base_enzyme_name, [sbo[count]])
                    current_prod_list  = [enz]
                    current_rxn_name = str(3*count)
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'binding'])
                    if debug:
                        print (current_rxn_name, current_react_list, current_prod_list)

                    #catalytic
                    current_react_list = [enz]
                    enz = enzyme_complex_name(current_base_enzyme_name, [pro[count]])
                    current_prod_list = [enz]
                    current_rxn_name = str(3*count+1)
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'catalytic'])
                    if debug:
                        print (current_rxn_name, current_react_list, current_prod_list)

                    # product release

                    current_react_list = [enz]
                    if (count+1) < len(sbo):
                        enz = f"{base_enzyme_name}_{count+1}"
                    else:
                        enz = base_enzyme_name
                    # correct: we want the following current one 
                    # unless it is the final one
                    current_prod_list = [enz,pro[count]]
                    current_rxn_name = str(3*count+2)
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'release'])
                    if debug:
                        print (current_rxn_name, current_react_list, current_prod_list)
                    
                count_start = 3*len(sbo) # check
                
        elif mechanism.lower() == 'rnd':
        #if mechanism.lower() == 'seq':
            # WARNING: implmentation is incomplete temporary stub
            #
            # TODO: finish implementation

            # random has all permutations
            # note that constraints on the kinetic parameters are not set here; only the network is

            # e.g. A + B -> P + Q via E with sbo: A, B & pro: P, Q
            # _0 E + A -> EA : binding
            # _1 EA + B -> EAB : binding
            # _2 E + B -> EB : binding
            # _3 EB + A -> EAB : binding
            # _4 EAB -> EPQ : catalytic
            # _5 EPQ -> EQ + P : release
            # _6 EQ -> E + Q : relase
            # _7 EPQ -> EP + Q : release
            # _8 EP -> E + P : release

            # that is, sbo [(0, 1), (1, 0)] and pro [(0, 1), (1, 0)]

            from itertools import permutations

            # substrate binding
            if sbo:
                # note: the order adds to the right in the order of the original sbo

                # use list and sort to ensure defined reproducible order
                sbo_perm = list(permutations(sbo))
                sbo_perm.sort()

                # need to track which forms are already written

            if pro:
                # note: the order adds to the right in the order of the original sbo

                # use list and sort to ensure defined reproducible order
                pro_perm = list(permutations(pro))
                pro_perm.sort()
        
        
        # inhibition steps
        reg_cumulative_count = 0
        
        imechanism = rxn_df['CI']
        if debug:
            print (imechanism)
        if imechanism:
            if debug:
                print('competitive')
            imechanism = imechanism.split(';')
            # currently K-FIT models CI as only binding the bare enzyme, 
            #   even for cases in which there are multiple items in SBO or 
            #   multiple items in CI list.
            
            # TODO:
            #   Here, all are written. For comparisons with K-FIT, 
            #   all after the first should have k parameters fixed to 0.
            for count,item in enumerate(imechanism):
                item_i = f"{item}_ci"
                enz = base_enzyme_name
                current_react_list = [enz,item]
                enz = enzyme_complex_name(base_enzyme_name,[item_i])
                current_prod_list  = [enz]
                current_rxn_name = f"i{reg_cumulative_count}"
                
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'competitive'])
                reg_cumulative_count += 1
        
        imechanism = rxn_df['UCI']
        if imechanism:
            if debug:
                print('uncompetitive')
            # K-FIT only has binding to the enzyme complex with all 
            #   substrates bound.
            
            # TODO:
            #   Here, all are written. For comparisons with K-FIT, 
            #   all after the first should have k parameters fixed to 0.
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                
                item_i = f"{item}_ui"
                enz = enzyme_complex_name(base_enzyme_name,sbo) # check
                current_react_list = [enz,item]
                enz = enzyme_complex_name(enz,[item_i])
                current_prod_list  = [enz]
                current_rxn_name = f"i{reg_cumulative_count}"
                
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'uncompetitive'])
                reg_cumulative_count += 1
                
        imechanism = rxn_df['NCI']
        if imechanism:
            if debug:
                print('noncompetitive')
            # need to check if need to add reaction between EI -> ESI at same rate as
                #    E + S -> ES
                # for now, only have E+I <=> EI and ES + I <=> ESI
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                item_i = f"{item}_ci"
                enz = base_enzyme_name
                current_react_list = [enz,item]
                enz = enzyme_complex_name(base_enzyme_name,[item_i])
                current_prod_list  = [enz]
                current_rxn_name = "i"+str(reg_cumulative_count)
                
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'noncompetitive'])
                reg_cumulative_count += 1
                
                item_i = f"{item}_ui"
                enz = enzyme_complex_name(base_enzyme_name,sbo) # check
                current_react_list = [enz,item]
                enz = enzyme_complex_name(enz,[item_i])
                current_prod_list  = [enz]
                current_rxn_name = f"i{reg_cumulative_count}"
                
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'elemental', 'noncompetitive'])
                reg_cumulative_count += 1
            
        # activation steps
        # for activators, in K-FIT the reaction is formed as
        #   A + E* <=> E
        #   where E* is the inactive enzyme and E is the active enzyme
        #   
        #   It is not clear that a system can handle a 'cold start', if a reaction
        #     requires a downstream metabolite to activate a reaction.
        #   presumably that won't be an issue
        #
        #   note: currently can only handle one activator per reaction
        imechanism = rxn_df['act']
        if imechanism:
            if debug:
                print('activation')
            imechanism = imechanism.split(';')
        # TODO : finish
        
        # generate dataframe
        df_new_reactions = df_join_reactions(pd.DataFrame(), rxn_list)
        
        return df_new_reactions

    def generate_mm_reactions(rxn_df, base_enzyme_name=None):
        mechanism = rxn_df.mechanism
        base_rxn_name = rxn_df.ID
        
        # TODO: validate all metabolites included

        rxn_list = []

        df_new_reactions = df_join_reactions()
  
        # substrates
        sbo = rxn_df['SBO']
        if sbo:
            sbo = sbo.split(';')
        # product
        pro = rxn_df['PRO']
        if pro:
            pro = pro.split(';')
        if debug:
            print (sbo,pro)

        if mechanism.lower() == 'seq':
            # substrate binding
            if sbo:
                # note: the order adds to the right
                for count,item in enumerate(sbo):
                    current_react_list = [item]
                    current_prod_list  = None
                    current_rxn_name = str(count)
                    
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'michaelis-menten', 'reactant'])
                    #print (current_rxn_name, current_react_list, current_prod_list)

            # product release
            if pro:
                # note: the order on the right pops off from left, which is the
                #   opposite of the add to the right for binding
                count_start = len(sbo)
                for count,item in enumerate(pro):
                    current_react_list = None 
                    current_prod_list = [item]
                    current_rxn_name = str(count+count_start)
                    
                    # add new data frame row
                    rxn_list.append([rxn_df.ID, current_rxn_name,
                                     current_react_list, current_prod_list,
                                     'michaelis-menten', 'product'])
        
        # inhibition steps
        reg_cumulative_count = 0

        imechanism = rxn_df['CI']
        if debug:
            print (imechanism)
        if imechanism:
            if debug:
                print('competitive')
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                current_react_list = [item]
                current_prod_list  = None
                current_rxn_name = f"i{reg_cumulative_count}"
                   
                # add new data frame row
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                current_react_list, current_prod_list,
                                'michaelis-menten', 'competitive'])
                reg_cumulative_count += 1

        imechanism = rxn_df['UCI']
        if debug:
            print (imechanism)
        if imechanism:
            if debug:
                print('uncompetitive')
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                current_react_list = [item]
                current_prod_list  = None
                current_rxn_name = f"i{reg_cumulative_count}"
                   
                # add new data frame row
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                current_react_list, current_prod_list,
                                'michaelis-menten', 'umcompetitive'])
                reg_cumulative_count += 1


        imechanism = rxn_df['NCI']
        if debug:
            print (imechanism)
        if imechanism:
            if debug:
                print('noncompetitive')
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                current_react_list = [item]
                current_prod_list  = None
                current_rxn_name = "i"+str(reg_cumulative_count)
                   
                # add new data frame row
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                current_react_list, current_prod_list,
                                'michaelis-menten', 'noncompetitive'])
                reg_cumulative_count += 1

        imechanism = rxn_df['act']
        if debug:
            print (imechanism)
        if imechanism:
            if debug:
                print('activator')
            imechanism = imechanism.split(';')
            for count,item in enumerate(imechanism):
                current_react_list = [item]
                current_prod_list  = None
                current_rxn_name = f"i{reg_cumulative_count}"
                   
                # add new data frame row
                rxn_list.append([rxn_df.ID, current_rxn_name,
                                current_react_list, current_prod_list,
                                'michaelis-menten', 'activator'])
                reg_cumulative_count += 1

        # generate dataframe
        df_new_reactions = df_join_reactions(pd.DataFrame(), rxn_list)
        
        return df_new_reactions


    df_expanded_reactions = df_join_reactions()
    df_reactions = pd.DataFrame()
    

    # main loop of function that loops over the reactions
    if testing_flag:
        df_reactions = generate_elemental_reactions(df_mechanism.loc[50])
    else:
        for rxn in df_mechanism.index:
            if debug:
                print (df_mechanism.loc[rxn].tolist())
            if mech_type.lower() == 'elemental':
                _tmp_er = generate_elemental_reactions(df_mechanism.loc[rxn])
                #print(ww)
            elif mech_type.lower() == "mm" or mech_type.lower() == 'michaelis-menten':
                
                _tmp_er = generate_mm_reactions(df_mechanism.loc[rxn])
            else:
                print (f"Unknown mechanism type {mech_type}")
                _tmp_er = generate_mm_reactions(df_mechanism.loc[rxn])
            df_reactions = pd.concat([df_reactions, _tmp_er],ignore_index=True)
    
    return df_reactions


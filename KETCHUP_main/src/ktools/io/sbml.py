""" 
sbml functions

collects basic functions that output results as sbml items

""" 
import libsbml
import pyomo.core.base
import cobra.core
from typing import Union


def sbml_check(value: Union[None,int], message: str) -> None:
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    # adapted from https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/libsbml-python-creating-model.html
    if value == None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                     + 'LibSBML returned error code ' + str(value) + ': "' \
                     + libsbml.OperationReturnValue_toString(value).strip() + '"'
            raise SystemExit(err_msg)
    else:
        return

def check_reaction_compartment(m_model:cobra.core.model.Model, rxn_id: str) -> list:
    """Returns list of compartments for the metabolites in selected reaction
    """
    rxn = m_model.reactions.get_by_id(rxn_id)
    compartments = list(dict.fromkeys( 
    [m.compartment for m in m_model.reactions.get_by_id(rxn_id).metabolites] 
    ))
    return compartments

def make_valid_sid(sid:str) -> str:
    """Returns string that contains only characters valid for use in id attributes.
    Specifically,
    letter ::= 'a'..'z','A'..'Z'
    digit  ::= '0'..'9'
    idChar ::= letter | digit | '_'
    """
    # This code substitutes some items with predefined substituions
    
    # TODO: use re to substitue other characters to the ASCII or unicode value
    
    sid = sid.replace('+','__plus__')
    sid = sid.replace('[','__lbrack__')
    sid = sid.replace(']','__rbrack__')
    sid = sid.replace(':','__colon__')
    sid = sid.replace("'",'__prime__')
    
    return sid

def v_rate_to_sbml(rate:str, block_id:str) -> str:
    """Returns rate law converted into valid sbml id attributes for the variables
    """
    rate = rate[rate.find('  ==  ')+6:] # slice off left side
    rate = rate.replace(block_id+'.','') # remove block id
    
    # assume that each term is multipled with k being the first
    # TODO make general to allow for more complex constraints
    if rate:
        tmp = rate.split('*')

        # TODO make prefix items defined
        for i,item in enumerate(tmp):
            if item[0] == 'k':
                tmp[i] = f"K_{tmp[0][:2]}_{make_valid_sid(tmp[0][3:-1])}".strip()
            else:
                if item[0] == 'c':
                    tmp[i] = f"S_{make_valid_sid(tmp[i][2:-1])}".strip()
                elif item[0] == 'e':
                    tmp[i] = f"E_{make_valid_sid(tmp[i][2:-1])}".strip()
        rate = '*'.join(tmp)
    return rate

def create_sbml_kinetic_model(pyomo_model:pyomo.core.base.PyomoModel.ConcreteModel,
                              block_id:int = 1) -> str:
    """Returns complete SBML Level 3 model of the kinetic parameterization
    at the solution
    """
    # First create empty SBMLDocument object
 
    try:
        document = libsbml.SBMLDocument(3, 2)
        sbml_check(document.setSBOTerm(675),           'set document SBO term')
    except ValueError:
        raise SystemExit('Could not create SBMLDocumention object')
    
    # Create basic Model object inside the SBMLDocument object
 
    model = document.createModel()
    sbml_check(model,                              'create model')
    sbml_check(model.setName(f"{pyomo_model.name}"),'set model name')
    sbml_check(model.setTimeUnits('second'),       'set model-wide time units')
    sbml_check(model.setExtentUnits('mole'),       'set model units of extent')
    sbml_check(model.setSubstanceUnits('mole'),    'set model substance units')
    
    # Create unit definitions

    per_second = model.createUnitDefinition()
    sbml_check(per_second,                         'create unit definition')
    sbml_check(per_second.setIdAttribute('per_second'),     'set unit definition id')
    unit = per_second.createUnit()
    sbml_check(unit,                               'create unit on per_second')
    sbml_check(unit.setKind(libsbml.UNIT_KIND_SECOND),     'set unit kind')
    sbml_check(unit.setExponent(-1),               'set unit exponent')
    sbml_check(unit.setScale(0),                   'set unit scale')
    sbml_check(unit.setMultiplier(1),              'set unit multiplier')

    # Create compartment inside this model

    for comp in pyomo_model.m_model.compartments:
        _tmp_c = model.createCompartment()
        sbml_check(_tmp_c,                                 'create compartment')
        sbml_check(_tmp_c.setIdAttribute(f"{comp}"),                'set compartment id')
        sbml_check(_tmp_c.setConstant(True),               'set compartment "constant"')
        sbml_check(_tmp_c.setSize(1),                      'set compartment "size"')
        sbml_check(_tmp_c.setSpatialDimensions(3),         'set compartment dimensions')
        sbml_check(_tmp_c.setUnits('litre'),               'set compartment size units')

    # Output WT solution
    # TODO: extend to other blocks
    
    cur_exp = pyomo_model.block_list.at(1) # pyomo_model.experiment0
    
    # Create species inside the model
    
    # First compounds
    
    for metab in pyomo_model.m_model.metabolites:
        _tmp_s = model.createSpecies()
        sbml_check(_tmp_s,                                 f'create species {metab.id}')
        sbml_check(_tmp_s.setIdAttribute(f"S_{make_valid_sid(metab.id)}"),          'set species id')
        sbml_check(_tmp_s.setName(f"{metab.name}"),          'set species name')
        sbml_check(_tmp_s.setCompartment(f"{metab.compartment}"), 'set species compartment')
        sbml_check(_tmp_s.setConstant(False),              'set "constant" attribute on species')
        sbml_check(_tmp_s.setInitialAmount(cur_exp.c[metab.id].value), 'set initial amount for species')
        sbml_check(_tmp_s.setSubstanceUnits('mole'),       'set substance units for species')
        sbml_check(_tmp_s.setBoundaryCondition(False),     'set "boundaryCondition" on species')
        sbml_check(_tmp_s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on species')
        sbml_check(_tmp_s.setSBOTerm(247),           'set species SBO term')
    # Next enzymes and complexes
    
    # Default compartment, using first defined in model
    # TODO: revisit, with possible list of comparment hierarchy
    default_compartment = list(pyomo_model.m_model.compartments)[0]
    
    for enz in cur_exp.e:
        _tmp_s = model.createSpecies()
        sbml_check(_tmp_s,                                 f'create species {enz}')
        sbml_check(_tmp_s.setIdAttribute(f"E_{make_valid_sid(enz)}"),          'set species id')
        sbml_check(_tmp_s.setName(f"{enz}"),          'set species name')
        sbml_check(_tmp_s.setCompartment(default_compartment), 'set species compartment')
        sbml_check(_tmp_s.setConstant(False),              'set "constant" attribute on species')
        sbml_check(_tmp_s.setInitialAmount(cur_exp.e[enz].value), 'set initial amount for species')
        sbml_check(_tmp_s.setSubstanceUnits('mole'),       'set substance units for species')
        sbml_check(_tmp_s.setBoundaryCondition(False),     'set "boundaryCondition" on species')
        sbml_check(_tmp_s.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on species')
        
        if '__plus__' not in _tmp_s.getIdAttribute():
            sbml_check(_tmp_s.setSBOTerm(14),           'set species SBO term')
        else:
            sbml_check(_tmp_s.setSBOTerm(253),           'set species SBO term')

    # Create parameter objects inside the model

    # TODO: Make general to allow all supported mechanism types
    
    for item in pyomo_model.kf:
        _tmp_k = model.createParameter()
        sbml_check(_tmp_k,                                  'create parameter k')
        sbml_check(_tmp_k.setIdAttribute(f"K_kf_{make_valid_sid(item)}"),            'set parameter k id')
        sbml_check(_tmp_k.setConstant(True),                'set parameter k "constant"')
        sbml_check(_tmp_k.setValue(pyomo_model.kf[item].value),'set parameter k value')
        sbml_check(_tmp_k.setUnits('per_second'),           'set parameter k units')
        sbml_check(_tmp_k.setSBOTerm(9),           'set parameter k SBO term')
        
    for item in pyomo_model.kr:
        _tmp_k = model.createParameter()
        sbml_check(_tmp_k,                                  'create parameter k')
        sbml_check(_tmp_k.setIdAttribute(f"K_kr_{make_valid_sid(item)}"),            'set parameter k id')
        sbml_check(_tmp_k.setConstant(True),                'set parameter k "constant"')
        sbml_check(_tmp_k.setValue(pyomo_model.kr[item].value),'set parameter k value')
        sbml_check(_tmp_k.setUnits('per_second'),           'set parameter k units')
        sbml_check(_tmp_k.setSBOTerm(9),           'set parameter k SBO term')
    
    # TODO: make SBO terms for parameters more specific
    
    # Create reactions inside the model
    
    for index, row in pyomo_model.mech_df.iterrows():
        _tmp_r = model.createReaction()
        sbml_check(_tmp_r,                                 'create reaction')
        sbml_check(_tmp_r.setIdAttribute(f"R_{make_valid_sid(row['rxn ID'])}_{make_valid_sid(row['step ID'])}"), 'set reaction id')
        sbml_check(_tmp_r.setName(pyomo_model.m_model.reactions.get_by_id(row['rxn ID']).name), 'set reaction name')
        sbml_check(_tmp_r.setReversible(True),            'set reaction reversibility flag')

        if row['mechanism'] == 'elemental':
            if row['type'] == 'binding':
                sbml_check(_tmp_r.setSBOTerm(177),           'set reaction SBO term')
            elif row['type'] == 'release':
                sbml_check(_tmp_r.setSBOTerm(180),           'set reaction SBO term')
            elif row['type'] == 'catalytic':
                sbml_check(_tmp_r.setSBOTerm(182),           'set reaction SBO term')
                # TODO: subdivide into those that are translocation (185), exchange (627)
                #       or remain as catalytic because of covalent bond changes
            else:
                sbml_check(_tmp_r.setSBOTerm(176),           'set reaction SBO term')
        else:
            sbml_check(_tmp_r.setSBOTerm(176),           'set reaction SBO term')

        # Set the reactants and products
        
        # TODO: revisit to combine into function
        
        for species in row['reactant']:
            _tmp_species = _tmp_r.createReactant()
            sbml_check(_tmp_species,                       'create reactant')
            #print (species)
            if species in pyomo_model.m_model.metabolites:
                sbml_check(_tmp_species.setSpecies(f"S_{make_valid_sid(species)}"),      'assign reactant species')
            else:
                sbml_check(_tmp_species.setSpecies(f"E_{make_valid_sid(species)}"),      'assign reactant species')  
                # Update enzyme compartment if all metabolites in the same compartment
                _rxn_compartments = check_reaction_compartment(pyomo_model.m_model,row['rxn ID'])
                if len(_rxn_compartments) == 1:
                    _tmp_s = model.getSpecies(f"E_{make_valid_sid(species)}")
                    sbml_check(_tmp_s.setCompartment(_rxn_compartments[0]), 'set species compartment')
            sbml_check(_tmp_species.setConstant(True),     'set "constant" on species')
        for species in row['product']:
            _tmp_species = _tmp_r.createProduct()
            sbml_check(_tmp_species,                       'create reactant')
            #print (species)
            if species in pyomo_model.m_model.metabolites:
                sbml_check(_tmp_species.setSpecies(f"S_{make_valid_sid(species)}"),      'assign product species')
            else:
                sbml_check(_tmp_species.setSpecies(f"E_{make_valid_sid(species)}"),      'assign product species')  
                # Update enzyme compartment if all metabolites in the same compartment
                _rxn_compartments = check_reaction_compartment(pyomo_model.m_model,row['rxn ID'])
                if len(_rxn_compartments) == 1:
                    _tmp_s = model.getSpecies(f"E_{make_valid_sid(species)}")
                    sbml_check(_tmp_s.setCompartment(_rxn_compartments[0]), 'set species compartment')
            sbml_check(_tmp_species.setConstant(True),     'set "constant" on species')

            # Set the reaction rate expression (the SBML "kinetic law")
            rxn_id = f"{row['rxn ID']}_{row['step ID']}"
            
            output_vf = v_rate_to_sbml(str(pyomo_model.vf_rate[cur_exp,rxn_id].expr), cur_exp.name)
            if output_vf == '':
                output_vf = '0'
        
            output_vr = v_rate_to_sbml(str(pyomo_model.vr_rate[cur_exp,rxn_id].expr), cur_exp.name)
            if output_vr == '':
                output_vr = '0'
        
            output_vnet = output_vf.strip() + ' - ' + output_vr.strip()

            _tmp_math_ast = libsbml.parseL3Formula(output_vnet)
            sbml_check(_tmp_math_ast,                      'create AST for rate expression')

            kinetic_law = _tmp_r.createKineticLaw()
            sbml_check(kinetic_law,                        'create kinetic law')
            sbml_check(kinetic_law.setMath(_tmp_math_ast), 'set math on kinetic law')
            sbml_check(kinetic_law.setSBOTerm(12),           'set kinetic law SBO term')
    
    return libsbml.writeSBMLToString(document)



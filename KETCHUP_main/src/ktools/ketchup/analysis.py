"""
analysis function

collects basic analysis functions for KETCHUP models

"""
import pandas as pd
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
import os
from typing import Union
import pyomo.core.base

def infeasible_constraints(model: pyomo.core.base.PyomoModel.ConcreteModel,
                           filename: str, tol: float = 10**-3):
    """
    Identifies and reports infeasible constraints in a provided Pyomo model.

    Parameters
    ----------
    model : pyomo.core.base.PyomoModel.ConcreteModel
        The Pyomo model to analyze.
    filename : str
        The path to the file used to append output reports.
    tol : float, optional
        The tolerance value used when considering a constraint as infeasible.
        Defaults to 1E-3.
    Returns
    -------
    None
    """

    with open(f"{filename}",'a+') as f:
        violations = []
        constraint_list = []
        for constr in model.component_data_objects(ctype=Constraint,active=True):
          try:
              c_val = value(constr)
              if c_val < constr.lb: 
                  if abs(c_val - constr.lb) > tol:
                      violation_str = f"lb violation - {constr}, {c_val}, {constr.lb}"
                      print(violation_str)              
                      f.write(f"{violation_str}\n")
                  violations.append( abs(c_val - constr.lb) )
                  constraint_list.append(constr)
              if c_val > constr.ub:
                  if abs(c_val - constr.ub) > tol:
                      violation_str = f"ub violation - {constr}, {c_val}, {constr.ub}"
                      print(violation_str) 
                      f.write(f"{violation_str}\n")
                  violations.append( abs(c_val - constr.ub) )
                  constraint_list.append(constr)
                  
          except ValueError:
              f.write(f"value error , {constr}\n")
        violation_str = f"max_violation - {max(violations)} | constraint {constraint_list[violations.index(max(violations))]}"
        print(violation_str)
        f.write(f"{violation_str}\n")
    return

def evaluate_stability(mech_df: pd.DataFrame, experiments: list, solnum: int,
                       pyomo_model: pyomo.core.base.PyomoModel.ConcreteModel,
                       time_solved: str, status: str) -> None:
    """
    Evaluates the stability of a parameterization of a kinetic model based on Jacobian analysis.

    Parameters
    ----------
    mech_df : pd.DataFrame
        DataFrame containing reaction and enzyme mechanism information.
    experiments : list
        List of experiment data blocks.
    solnum : int
        Solution number.
    pyomo_model : pyomo.core.base.PyomoModel.ConcreteModel
        The Pyomo model to analyze.
    time_solved : str
        Time taken to solve the model. Used to create unique filename of output.
    status : str
        Solver status of the solution.

    Returns
    -------
    None
    """

    nlp = PyomoNLP(pyomo_model)
    m = nlp.evaluate_jacobian().toarray() 
    #h = nlp.evaluate_hessian_lag().todense()
    #g = nlp.evaluate_grad_objective()

    elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
    rki = {}
    for i,exp in enumerate(experiments):
        for j,es in enumerate(pyomo_model.rxn_enz_sum.keys()):
            rki.update({(i*len(pyomo_model.rxn_enz_sum.keys())+j):es})
            rki.update({len(experiments)*len(pyomo_model.rxn_enz_sum.keys())+(i*len(pyomo_model.rxn_enz_sum.keys())+j):es})
    vf_ind = nlp.get_constraint_indices([pyomo_model.vf_rate])
    vr_ind = nlp.get_constraint_indices([pyomo_model.vr_rate])
    rxn_count = len(pyomo_model.ELEMENTALSTEP_F)
    stability = 'stable'
    for e,exp in enumerate(experiments):
        J = m[vf_ind[e*rxn_count:(e+1)*rxn_count]+vr_ind[e*rxn_count:(e+1)*rxn_count],:][:,nlp.get_primal_indices([pyomo_model.kf])+nlp.get_primal_indices([pyomo_model.kr])]
    
        zero_eig = []
        thres = 10**-3 # threshold
        w = np.linalg.eigvals(J)
        p = 0;    n = 0;    z = 0
        for j,i in enumerate(w):
            if i > thres  : 
               p += 1
               zero_eig.append((j,rki[j]))
          
            if i == 0: z += 1
            if i < -thres  : n += 1
        if p > 0: stability = 'unstable'
        #print(f"{exp} - pos {p}\nneg {n}\nzer {z}",flush=True)    

    if status == 'optimal':
        result_dump(f"s_{status}_{stability}_results",solnum,pyomo_model,time_solved,status)
    return None



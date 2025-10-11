"""
outputs function

collects basic functions that output results

"""
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
import os
from typing import Union


def result_dump(name: str, solnum: int, model, time_solved: float, status: str,
                data_type_category: str = 'static',
                mechanism_type_category: str = 'elemental') -> None:
    """
    Writes solution results to a text file in JSON format.

    Parameters
    ----------
    name : str
        The base name of the desired result output file.
    solnum : int
        The solution identifier number.
    model : pyomo.opt.results.results_.SolverResults
        The Pyomo ConcreteModel representing the model.  Contains the data to be written.
    time_solved : float
        The time taken to solve the model.
    status : str
        The solver status (e.g., 'optimal').
    data_type_category : str
        Category of data ('static' or 'dynamic').
    mechanism_type_category : str
        Category of rate law mechanism.

    Returns
    -------
    None
    """

    import json
    
    fn = f"{name}_{str(solnum)}.txt"
    res = {} # results

    if mechanism_type_category == 'elemental':
        res_kf = {}
        for item in model.kf:
            res_kf[item] = model.kf[item].value
        res['kf'] = res_kf # forward k

        res_kr = {}
        for item in model.kr:
            res_kr[item] = model.kr[item].value
        res['kr'] = res_kr # reverse k
    elif mechanism_type_category == 'custom':
        res_kp = {}
        for item in model.kinetic_parameters:
            _tmp_kp = {}
            for k in model.kinetic_parameters[item]:
                _tmp_kp.update({k:model.kinetic_parameters[item][k].value})
            res_kp.update({item:_tmp_kp})
        res['kp'] = res_kp # k_p
    
    exp_keys = list(model.data.keys())
    
    res_c = {} # compounds (metabolites)
    res_e = {} # enzymes
    res_rate = {} # kinetic rates
    res_SSR = {} # sum of square residuals
    res_time = {} # time
    total_error = 0

    # loop over experiment blocks
    for count, key in enumerate(exp_keys):
        #try:
        _tmp_c = {}
        _tmp_e = {}
        _tmp_rate = {}

        if count == 0:
            cur_exp = model.experiment0
        else:
            cur_exp = eval(f"model.experiment{str(count)}")

        if data_type_category == 'static':
            for item in cur_exp.c:
                _tmp_c[str(item)] = cur_exp.c[item].value
            for item in cur_exp.e:
                _tmp_e[str(item)] = cur_exp.e[item].value
            for item in cur_exp.rate:
                _tmp_rate[str(item)] = cur_exp.rate[item].value
            total_error += cur_exp.error.value
        elif data_type_category == 'dynamic':
            for item in cur_exp.c:
                _tmp_c[str(item)] = cur_exp.c[item].value
            for item in cur_exp.e:
                _tmp_e[str(item)] = cur_exp.e[item].value
            for item in cur_exp.rate:
                _tmp_rate[str(item)] = cur_exp.rate[item].value
            total_error += cur_exp.error.value
        else:
            for item in cur_exp[10].c:
                _tmp_c[str(item)] = cur_exp[10].c[item].value
            for item in cur_exp[10].e:
                _tmp_e[str(item)] = cur_exp[10].e[item].value
            for item in cur_exp[10].rate:
                _tmp_rate[str(item)] = cur_exp[10].rate[item].value
            total_error += cur_exp.error.value
        res_c[key] = _tmp_c
        res_e[key] = _tmp_e
        res_rate[key] = _tmp_rate
    res['c'] = res_c
    res['e'] = res_e
    res['rate'] = res_rate
    res['SSR'] = total_error 
    res['time'] = time_solved

    with open(fn, 'w') as json_file:
        # json.dump(res, json_file)
        json.dump(res, json_file, sort_keys=True, indent='\t')
        print(f"Successful export of data into {fn}")

    if status == 'optimal':
        with open('optimal_solutions.txt', 'a+') as f:
            f.write(f"iteration {solnum} - SSR {round(float(total_error),3)} - obj value - {round(float(model.obj()),3)} time - {time_solved}\n")

    with open ('total_runs.txt','a+') as f:
        f.write(f"iteration:{solnum} - status:{status} - time:{time_solved}\n")   
    return None


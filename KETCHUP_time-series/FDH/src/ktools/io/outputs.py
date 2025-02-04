"""
outputs function

collects basic functions that output results

"""
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
import os
from typing import Union


def result_dump(name, solnum, model, time_solved, status, static=1) -> None:
    import json
    
    fn = f"{name}_{str(solnum)}.txt"
    res = {}
    m_type = model.options['mechanism_type']
    if m_type == 'elemental':
        res_kf = {}
        for item in model.kf:
            res_kf[item] = model.kf[item].value
        res['kf'] = res_kf
        
        res_kr = {}
        for item in model.kr:
            res_kr[item] = model.kr[item].value
        res['kr'] = res_kr
    elif m_type == 'custom':
        res_kp = {}
        for item in model.kinetic_parameters:
            _tmp_kp = {}
            for k in model.kinetic_parameters[item]:
                _tmp_kp.update({k:model.kinetic_parameters[item][k].value})
            res_kp.update({item:_tmp_kp})
        res['kp'] = res_kp

    exp_keys = list(model.data.keys())
    
    res_c = {}
    res_e = {}
    res_rate = {}
    res_SSR = {}
    res_time = {}
    total_error = 0
    for count, key in enumerate(exp_keys):
        #try:
         _tmp_c = {}
         _tmp_e = {}
         _tmp_rate = {}

         if count == 0:
             cur_exp = model.experiment0
         else:
             cur_exp = eval(f"model.experiment{str(count)}")

         if static:
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
         res_c[key] = _tmp_c
         res_e[key] = _tmp_e
         res_rate[key] = _tmp_rate
    res['c'] = res_c
    res['e'] = res_e
    res['rate'] = res_rate
    res['SSR'] = total_error 
    res['time'] = time_solved

    with open(fn, 'w') as json_file:
        json.dump(res, json_file,sort_keys=True, indent='\t', separators=(',', ': '))
        print(f"sucessful export of data into {fn}")

    if status == 'optimal':
        with open('optimal_solutions.txt', 'a+') as f:
            f.write(f"iteration {solnum} - SSR {round(float(total_error),3)} - obj value - {round(float(model.obj()),3)} time - {time_solved}\n")

    with open ('totalruns.txt','a+') as f:
        f.write(f"iteration:{solnum} - status:{status} - time:{time_solved}\n")   
    return None



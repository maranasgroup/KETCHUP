"""
outputs function

print results as json

"""
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
import os


def result_dump(name,solnum,model,time_solved,status,static=1,):
    import json
    
    fn = name+'_'+str(solnum)+'.txt'
    res = {}
    
    res_kf = {}
    for item in model.kf:
        res_kf[item] = model.kf[item].value
    res['kf'] = res_kf
    
    res_kr = {}
    for item in model.kr:
        res_kr[item] = model.kr[item].value
    res['kr'] = res_kr
    
    exp_keys = list(model.data.keys())
    
    res_c = {}
    res_e = {}
    res_rate = {}
    res_SSR = {}
    res_time = {}
    total_error = 0
    for count,key in enumerate(exp_keys):
        #try:
         _tmp_c = {}
         _tmp_e = {}
         _tmp_rate = {}
         if count == 0:
             cur_exp = model.experiment0
         else:
             cur_exp = eval('model.experiment'+str(count))
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
        #except:
        #    continue
    res['c'] = res_c
    res['e'] = res_e
    res['rate'] = res_rate
    res['SSR'] = total_error 
    res['time'] = time_solved

    with open(fn, 'w') as json_file:
        json.dump(res, json_file)

    print(f"sucessful export of data into {fn}")


    if status == 'optimal':
        if os.path.isfile('optimal_solutions.txt'):
            f = open('optimal_solutions.txt','a')
        else:
            f = open('optimal_solutions.txt','w')
        f.write(f'iteration {solnum} - SSR {round(float(total_error),3)} - obj value - {round(float(model.obj()),3)} time - {time_solved}\n')

    if os.path.isfile('totalruns.txt'):
        f = open('totalruns.txt','a')
    else:
       f = open('totalruns.txt','w')
    f.write(f'iteration:{solnum} - status:{status} - time:{time_solved}\n')   
    return None

def infeasible_constraints(model,filename,tol=10**-3):
    with open(f'{filename}','a') as f:
        violations = []
        constraint_list = []
        for constr in model.component_data_objects(ctype=Constraint,active=True):
          try:
              c_val = value(constr)
              if c_val < constr.lb: 
                  if abs(c_val - constr.lb) > tol: 
                      f.write(f'lb violation - {constr}, {c_val}, {constr.lb}\n')
                      print(f'lb violation - {constr}, {c_val}, {constr.lb}')              
                  violations.append( abs(c_val - constr.lb) )
                  constraint_list.append(constr)
              if c_val > constr.ub:
                  if abs(c_val - constr.ub) > tol:
                      f.write(f'ub violation - {constr}, {c_val}, {constr.ub}\n')
                      print(f'ub violation - {constr}, {c_val}, {constr.ub}') 
                  violations.append( abs(c_val - constr.ub) )
                  constraint_list.append(constr)
                  
          except ValueError:
              f.write(f'value error , {constr}\n')
        print(f'max_violation - {max(violations)} | constraint {constraint_list[violations.index(max(violations))]}')
        f.write(f'max_violation - {max(violations)} | constraint {constraint_list[violations.index(max(violations))]}\n')
        f.close()



def evaluate_stability(mech_df,experiments,solnum,dae_model,time_solved,status):
    nlp = PyomoNLP(dae_model)
    m = nlp.evaluate_jacobian().toarray() 
    #h = nlp.evaluate_hessian_lag().todense()
    #g = nlp.evaluate_grad_objective()


    elemental_steps = [f"{r}_{mech_df['step ID'].values[i]}" for i,r in enumerate(mech_df['rxn ID'].tolist())]
    rki = {}
    for i,exp in enumerate(experiments):
        for j,es in enumerate(dae_model.rxn_enz_sum.keys()):
            rki.update({(i*len(dae_model.rxn_enz_sum.keys())+j):es})
            rki.update({len(experiments)*len(dae_model.rxn_enz_sum.keys())+(i*len(dae_model.rxn_enz_sum.keys())+j):es})
    vf_ind = nlp.get_constraint_indices([dae_model.vf_rate])
    vr_ind = nlp.get_constraint_indices([dae_model.vr_rate])
    rxn_count = len(dae_model.ELEMENTALSTEP_F)
    stability = "stable"
    for e,exp in enumerate(experiments):
        J = m[vf_ind[e*rxn_count:(e+1)*rxn_count]+vr_ind[e*rxn_count:(e+1)*rxn_count],:][:,nlp.get_primal_indices([dae_model.kf])+nlp.get_primal_indices([dae_model.kr])]
    
        zero_eig = []
        thres = 10**-3
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

    if status == "optimal": result_dump(f"s_{status}_{stability}_results",solnum,dae_model,time_solved,status)
    return None
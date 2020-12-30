

# codegen的模块套一层壳

import gfold_subproblem_solver
import numpy as np

import GFOLD_params

def wrap(v):
    return np.matrix([[v], [0]])


def solve(params, params_super = None, verbose=False):
    
    #super params
    if (params_super == None):
        params_super = GFOLD_params.SuperParams() # default
    N = params_super.N
    
    res = gfold_subproblem_solver.cg_solve(
        x0 = params.x0.reshape(6, 1),
        xf = params.xf.reshape(6, 1),
        z0_term_inv = params.z0_term_inv.reshape(1, N),
        z0_term_log = params.z0_term_log.reshape(1, N),
        g_vec = params.g.reshape(3, 1),
        p_cs_cos = params.p_cs_cos.reshape(1, N),
        m_wet_log = wrap(params.m_wet_log),
        sparse_params = np.array([
            params.alpha_dt, 
            params.G_max, 
            params.V_max, 
            params.y_gs_cot, 
            params.r1, 
            params.r2, 
            params.tf
            ]).reshape(7, 1)
    )
    #print(res[1]['status'])
    return (
        res,
        np.array(res[0]['var_x']), # r,v
        np.array(res[0]['var_u']), # u (acceleration)
        np.exp(np.array(res[0]['var_z'])) # mass
    ) if res[1]['status'] == 'optimal' else (None, None, None, None)


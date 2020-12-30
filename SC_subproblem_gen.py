
# codegen的模块套一层壳

import sc_subproblem_solver
import numpy as np

import SC_params

def wrap(v):
    return np.matrix([[v], [0]])

def solve(params, params_super = None):
    
    #super params
    if (params_super == None):
        params_super = SC_params.SuperParams() # default
    K = params_super.K
   
    
    res = sc_subproblem_solver.cg_solve(
        A=params.A.reshape(K * 14, 14),
        B=params.B.reshape(K * 14, 3),
        C=params.C.reshape(K * 14, 3),
        S=params.S.reshape(K * 14, 1),
        z=params.z.reshape(K * 14, 1),
        x_last=params.x_last,
        u_last=params.u_last,
        u_last_dir=params.u_last_dir,
        s_last=wrap(params.s_last),
        w_nu=wrap(params.w_nu),
        w_delta=wrap(params.w_delta),
        w_delta_s=wrap(params.w_delta_s),
        
        x_initial=params.x_initial,
        x_final=params.x_final,
        #sparse
        m_dry=wrap(params.m_dry),
        tan_gamma_gs=wrap(params.tan_gamma_gs),
        cos_theta_max=(np.array([params.cos_theta_max] * K) if type(params.cos_theta_max)!=np.ndarray else params.cos_theta_max).reshape(1, K),
        omega_max=wrap(params.omega_max),
        cos_delta_max=wrap(params.cos_delta_max),
        T_max=wrap(params.T_max),
        T_min=wrap(params.T_min)
    )
    return (res, 
        np.array(res[0]['x']), 
        np.array(res[0]['u']), 
        res[0]['s'][0,0], 
        np.array(res[0]['nu']), 
        np.array(res[0]['delta']), 
        res[0]['delta_s']) #tuple(结果，状态，控制，时间scale)

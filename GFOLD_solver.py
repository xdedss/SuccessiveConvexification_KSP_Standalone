

import numpy as np
import math
from time import time 

import GFOLD_params


def solve(vessel_profile, vessel_initial, vessel_final, params_super=None, use_c=False, verbose=False):
    #super params
    if (params_super == None):
        params_super = GFOLD_params.SuperParams() # default
    N = params_super.N
    
    # solver
    if (use_c):
        import GFOLD_problem_gen as solver
    else:
        import GFOLD_problem as solver
    
    params = GFOLD_params.Params(N)
    
    #initial and final
    params.x0[0:3, 0] = vessel_initial.pos
    params.x0[3:6, 0] = vessel_initial.vel
    params.xf[0:3, 0] = vessel_final.pos
    params.xf[3:6, 0] = vessel_final.vel
    
    dt = vessel_profile.time_guess / N
    alpha = 1 / 9.80665 / vessel_profile.isp
    alpha_dt = alpha * dt
    t = np.linspace(0, (N-1) * dt, N)
    params.r1 = vessel_profile.T_min
    params.r2 = vessel_profile.T_max
    z0_term = vessel_initial.mass - alpha * params.r2 * t
    params.z0_term_inv = (1 / z0_term).reshape(1, N)
    params.z0_term_log = np.log(z0_term).reshape(1, N)
    params.g = vessel_profile.g
    params.alpha_dt = alpha_dt
    params.G_max = 100 # todo
    params.V_max = 500 # todo
    params.y_gs_cot = 1 / np.tan(vessel_profile.gamma_gs)
    params.m_wet_log = np.log(vessel_initial.mass)
    params.tf = vessel_profile.time_guess
    
    if (type(vessel_profile.theta_max) == np.ndarray):
        n_space = np.linspace(0, N-1, len(vessel_profile.theta_max))
        p_cs_array = np.array([np.interp(n, n_space, vessel_profile.theta_max) for n in range(N)])
    else:
        p_cs_array = np.array([vessel_profile.theta_max] * N)
    params.p_cs_cos = np.cos(p_cs_array)
    
    
    
    obj_opt, x, u, m = solver.solve(params, verbose=verbose)
    return x, u, m

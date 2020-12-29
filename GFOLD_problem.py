# -*- coding: utf-8 -*-
# GFOLD_static_p3p4


min_=min
from cvxpy import *
import cvxpy_codegen as cpg
from time import time
import numpy as np
import sys

import GFOLD_params

''' As defined in the paper...

 PROBLEM 3: Minimum Landing Error (tf roughly solved)
 MINIMIZE : norm of landing error vector
 SUBJ TO  :
            0) initial conditions satisfied (position, velocity)
            1) final conditions satisfied (altitude, velocity)
            2) dynamics always satisfied
            3) x stays in cone at all times
            4) relaxed convexified mass and thrust constraints
            5) thrust pointing constraint
            6) sub-surface flight constraint

 PROBLEM 4: Minimum Fuel Use
 MAXIMIZE : landing mass, opt variables are dynamical and
 SUBJ TO  :
            0) same constraints as p1, plus:
            1) landing point must be equal or better than that found by p1

'''


def solve(params, params_super = None, codegen = False, verbose=False):
    #super params
    if (params_super == None):
        params_super = GFOLD_params.SuperParams() # default
    N = params_super.N
    
    #优化变量
    x =Variable(6,N,name='var_x') # state vector (3position,3velocity)
    u =Variable(3,N,name='var_u') # u = Tc/mass because Tc[:,n]/m[n] is not allowed by DCP
    z= Variable(1,N,name='var_z')  # z = ln(mass)
    s= Variable(1,N,name='var_s') # thrust slack parameter
    
    # Parameters
    x0 = Parameter(6, 1, name="x0")
    xf = Parameter(6, 1, name="xf")
    z0_term_inv = Parameter(1, N, name="z0_term_inv", sign='positive')
    z0_term_log = Parameter(1, N, name="z0_term_log")
    g = Parameter(3, 1, name="g_vec")
    p_cs_cos = Parameter(1, N, name='p_cs_cos')
    sparse_params = Parameter(7, 1, name="sparse_params", sign='positive')
    m_wet_log = Parameter(2, 1, name='m_wet_log')
    if (not codegen):
        x0.value = params.x0.reshape(6, 1)
        xf.value = params.xf.reshape(6, 1)
        z0_term_inv.value = params.z0_term_inv.reshape(1, N)
        z0_term_log.value = params.z0_term_log.reshape(1, N)
        g.value = params.g.reshape(3, 1)
        p_cs_cos.value = params.p_cs_cos.reshape(1, N)
        m_wet_log.value = [params.m_wet_log, 0]
        sparse_params.value = np.array([
            params.alpha_dt, 
            params.G_max, 
            params.V_max, 
            params.y_gs_cot, 
            params.r1, 
            params.r2, 
            params.tf
            ]).reshape(7, 1)
    alpha_dt, G_max, V_max, y_gs_cot, r1, r2, tf_ = sparse_params

    dt = tf_ * (1/N)  # Integration dt
    
    # constraints
    con = [] 
    
    con += [x[0:3,0]  == x0[0:3]] # initial pos
    con += [x[3:6,0]  == x0[3:6]] # initial vel
        
    con += [x[0:3,N-1] == xf[0:3]] # final pos
    con += [x[3:6,N-1]== xf[3:6]] # final vel

    con += [s[0,N-1] == 0] # thrust at the end must be zero
    con += [u[:,0] == s[0,0]*np.array([1,0,0])] # thrust direction starts straight
    con += [u[:,N-1] == s[0,N-1]*np.array([1,0,0])] # and ends straight
    con += [z[0,0] == m_wet_log[0,0]] # convexified (7)

    
    for n in range(0,N-1):
        
        #dynamics
        con += [x[3:6,n+1] == x[3:6,n] + (dt*0.5)*((u[:,n]+g[:,0]) + (u[:,n+1]+g[:,0]))]
        con += [x[0:3,n+1] == x[0:3,n] + (dt*0.5)*(x[3:6,n+1]+x[3:6,n])]

        # glideslope cone
        con += [ norm( (x[0:3,n])[1:3] ) - y_gs_cot*(x[0,n])  <= 0 ]
            
        con += [ norm(x[3:6,n]) <= V_max ] # velocity
        #con += [norm(u[:,n+1]-u[:,n]) <= dt*T_max/m_dry * 3]
        con += [z[0,n+1] == z[0,n] - (alpha_dt*0.5)*(s[0,n] + s[0,n+1])] # mass decreases
        con += [norm(u[:,n]) <= s[0,n]] # limit thrust magnitude & also therefore, mass

        # Thrust pointing constraint
        con += [ u[0,n] >= p_cs_cos[0,n]*s[0,n]  ]

        if n > 0:
            
            #z0_term = m_wet - alpha * r2 * (n) * dt  # see ref [2], eq 34,35,36
            
            #z0 = log(z0_term)
            
            z0 = z0_term_log[0,n]
            
            mu_1 = r1*(z0_term_inv[0,n])
            mu_2 = r2*(z0_term_inv[0,n])
            
            #更正一处原项目与论文不符之处
            # 示意图：https://www.desmos.com/calculator/wtcfgnepe1
            con += [s[0,n] >= mu_1 * (1 - (z[0,n] - z0) + (z[0,n] - z0)**2 *0.5)] # lower thrust bound
            con += [s[0,n] <= mu_2 * (1 - (z[0,n] - z0))] # upper thrust bound


    #Objective
    objective = Minimize(-z[0,N-1])
    problem=Problem(objective, con)
    
    if codegen:
        cpg.codegen(problem, codegen_path)
    else:
        obj_opt = problem.solve(solver=ECOS, verbose=verbose)
        return (
            obj_opt,
            np.array(x.value), # r,v
            np.array(u.value), # u (acceleration)
            np.exp(np.array(z.value)) # mass
        )

if __name__ == '__main__':
    if (len(sys.argv) > 2 and sys.argv[1] == 'codegen'):
        codegen_path = sys.argv[2]
        solve(None, None, True)
    else:
        print("invalid input")
        print(sys.argv)
    
    

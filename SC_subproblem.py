
#subproblem: x^ v^ sigma^ => next x y sigma
#对应原论文 Problem 2. : Convex Discrete-Time Fixed-Final-Time Problem

from cvxpy import *
import cvxpy_codegen as cpg
from time import time
import numpy as np
import sys

import SC_params

#codegen时params为None，super必须提供 K
#超参数贯穿整个codegen
def solve(params, params_super = None, codegen = False):
    #super params
    if (params_super == None):
        params_super = SC_params.SuperParams() # default
    K = params_super.K
    
    #优化变量
    x =         Variable(14, K, name='x')
    u =         Variable(3, K, name='u')
    s =         Variable(1, 1, name='s')
    nu =        Variable(14, K-1, name='nu')
    delta =     Variable(1, K, name='delta')
    delta_s =   Variable(1, 1, name='delta_s')
    
    #参数
#    A =         [Parameter(14, 14, name='A_%s'%i) for i in range(K)]
#    B =         [Parameter(14, 3, name='B_%s'%i) for i in range(K)]
#    C =         [Parameter(14, 3, name='C_%s'%i) for i in range(K)]
#    S =         [Parameter(14, 1, name='S_%s'%i) for i in range(K)]
#    z =         [Parameter(14, 1, name='z_%s'%i) for i in range(K)]
    #暴力数组会超python255参数限制。。。
    A =         Parameter(14 * K, 14, name='A')
    B =         Parameter(14 * K, 3, name='B')
    C =         Parameter(14 * K, 3, name='C')
    S =         Parameter(14 * K, 1, name='S')
    z =         Parameter(14 * K, 1, name='z')
    x_last      = Parameter(14, K, name='x_last')
    u_last      = Parameter(3, K, name='u_last')
    u_last_dir      = Parameter(3, K, name='u_last_dir')
    s_last      = Parameter(2, 1, name='s_last') #2行1列，第二列占位用，因为1行1列在codegen会出bug
    w_nu        = Parameter(2, 1, name='w_nu', sign="positive") #权重非负，保证dcp
    w_delta     = Parameter(2, 1, name='w_delta', sign="positive")
    w_delta_s   = Parameter(2, 1, name='w_delta_s', sign="positive")
    
    x_initial = Parameter(14, 1, name='x_initial')
    x_final   = Parameter(14, 1, name='x_final')
    #sparse
    m_dry = Parameter(2, 1, name='m_dry')
    tan_gamma_gs = Parameter(2, 1, name='tan_gamma_gs', sign="positive")
    cos_theta_max = Parameter(2, 1, name='cos_theta_max')
    omega_max = Parameter(2, 1, name='omega_max')
    cos_delta_max = Parameter(2, 1, name='cos_delta_max', sign="positive")
    T_max = Parameter(2, 1, name='T_max')
    T_min = Parameter(2, 1, name='T_min')
    
    if (not codegen): #填入实际参数
        A.value = params.A.reshape(K * 14, 14)
        B.value = params.B.reshape(K * 14, 3)
        C.value = params.C.reshape(K * 14, 3)
        S.value = params.S.reshape(K * 14, 1)
        z.value = params.z.reshape(K * 14, 1)
        x_last.value = params.x_last
        u_last.value = params.u_last
        u_last_dir.value = params.u_last_dir
        s_last.value = [params.s_last, 0]
        w_nu.value = [params.w_nu, 0]
        w_delta.value = [params.w_delta, 0]
        w_delta_s.value = [params.w_delta_s, 0]
        
        x_initial.value = params.x_initial
        x_final.value = params.x_final
        #sparse
        m_dry.value = [params.m_dry, 0]
        tan_gamma_gs.value = [params.tan_gamma_gs, 0]
        cos_theta_max.value = [params.cos_theta_max, 0]
        omega_max.value = [params.omega_max, 0]
        cos_delta_max.value = [params.cos_delta_max, 0]
        T_max.value = [params.T_max, 0]
        T_min.value = [params.T_min, 0]
    
    #限制条件
    cons = []
    
#（1）边界条件
    #初始
    cons += [
        x[0, 0]     == x_initial[0, 0], # mass
        x[1:4, 0]   == x_initial[1:4, 0], # position
        x[4:7, 0]   == x_initial[4:7, 0], # velocity
        x[7:11, 0]  == x_initial[7:11, 0], # quanternion 更改：初态姿态固定
        x[11:14, 0] == x_initial[11:14, 0], # angular vel
    ]
    #结束
    cons += [
        # x[0, K-1]     == x_final[0, 0], # mass
        x[1:4, K-1]   == x_final[1:4, 0], # position
        x[4:7, K-1]   == x_final[4:7, 0], # velocity
        #x[7:11, K-1]  == x_final[7:11, 0], # quanternion
        x[11:14, K-1] == x_final[11:14, 0], # angular vel
    ]
    #thrust最后朝下
    cons += [
        u[1, K-1] == 0,
        u[2, K-1] == 0,
    ]
    
#（2）动力学方程
    for k in range(K - 1):
        cons += [
            x[:, k+1] == A[14*k:14*(k+1),:]*x[:, k] + B[14*k:14*(k+1),:]*u[:, k] + C[14*k:14*(k+1),:]*u[:, k+1] + S[14*k:14*(k+1),:]*s + z[14*k:14*(k+1),:] + nu[:, k]
        ]
    
#（3）状态限制
    for k in range(K):
        cons += [
            x[0, k] >= m_dry[0,0], #燃料耗尽
            norm(x[2:4, k]) * tan_gamma_gs[0,0] <= x[1, k], #在锥内部
            cos_theta_max[0,0] <= 1 - 2 * sum_squares(x[9:11, k]), # 倾角
            norm(x[11:14, k]) <= omega_max[0,0], #角速度
        ]
    cons += [0 == x[9:11, K-1]] # 规定最终头朝上但是滚转轴随意  即四元数jk分量为0
    #cons += [0 == x[8, 0]]
    cons += [s >= 0]
    
#（4）输入量（推力）限制
    for k in range(K):
        cons += [
            #T_min[0,0] <= u_last_dir[:, k].T * u[:, k], #最小推力线性近似
            T_min[0,0] <= u_last_dir[0, k]*u[0, k] + u_last_dir[1, k]*u[1, k] + u_last_dir[2, k]*u[2, k], #点乘展开写
            norm(u[:, k]) <= T_max[0,0], #最大推力凸约束
            norm(u[:, k]) * cos_delta_max[0,0] <= u[0, k], #gimbal限制
        ]
    
#（5）trust region
    for k in range(K):
        dx = x[:, k] - x_last[:, k]
        du = u[:, k] - u_last[:, k]
        cons += [
            sum_squares(dx) + sum_squares(du) <= delta[:, k]
        ]
    cons += [norm(s - s_last[0,0], 1) <= delta_s]
    
# Objective:
    objective = Minimize(
        s + w_nu[0,0]      * norm(nu, 1)  # virtual control（1范数）
          + w_delta[0,0]   * norm(delta)   # trust region on dynamics（2范数）
          + w_delta_s[0,0] * norm(delta_s, 1)  # trust region on sigma（1范数）
    )
    
# Problem
    problem = Problem(objective, cons)
    
    #print('is DCP: %s' % problem.is_dcp()) #检查是否符合凸优化规则

# Solve or Codegen
    if (codegen):
        cpg.codegen(problem, codegen_path)
    else:
        obj_opt = problem.solve(solver=ECOS, verbose=False)
        return (obj_opt, 
            np.array(x.value), 
            np.array(u.value), 
            s.value, 
            np.array(nu.value), 
            np.array(delta.value), 
            delta_s.value) #tuple(结果，状态，控制，时间scale)
    

if __name__ == '__main__':
    if (len(sys.argv) > 2 and sys.argv[1] == 'codegen'):
        codegen_path = sys.argv[2]
        solve(None, None, True)
    else:
        print("invalid input")
        print(sys.argv)
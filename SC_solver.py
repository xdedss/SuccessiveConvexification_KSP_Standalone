

import numpy as np
import math
from time import time
import pickle, sys

import SC_params
from SC_integrator import Integrator


def solve(vessel_profile, vessel_initial, vessel_final, solver_options=None, params_super=None, use_c = False, verbose = False):
    # default params
    if (params_super == None):
        params_super = SC_params.SuperParams()
    if (solver_options == None):
        solver_options = SC_params.SolverOptions()
    
    #normalize
    Ut = vessel_profile.time_guess / 5.
    Ul = vessel_initial.pos[0] / 20.
    Um = vessel_initial.mass / 1.5
    vessel_profile = SC_params.normalize(vessel_profile, Ut, Ul, Um)
    vessel_initial = SC_params.normalize(vessel_initial, Ut, Ul, Um)
    vessel_final = SC_params.normalize(vessel_final, Ut, Ul, Um)
    
    # solver
    if (use_c):
        import SC_subproblem_gen as solver
    else:
        import SC_subproblem as solver
    import GFOLD_solver, GFOLD_params
    
    # super params
    K = params_super.K
    dt = 1 / (K - 1)
    iterations = solver_options.iterations
    if (verbose):
        print('K=%s, iterations=%s' % (K, iterations))
    
    # 参数的转换 将飞行器参数计算转化后填入求解器参数
    params = SC_params.Params(K)
    #sparse
    params.m_dry = vessel_profile.m_dry
    params.tan_gamma_gs = math.tan(vessel_profile.gamma_gs)
    params.cos_theta_max = np.cos(vessel_profile.theta_max)
    params.omega_max = vessel_profile.omega_max
    params.cos_delta_max = math.cos(vessel_profile.delta_max)
    params.T_min = vessel_profile.T_min
    params.T_max = vessel_profile.T_max
    
    #initial
    params.x_initial[0:1, 0] = vessel_initial.mass
    params.x_initial[1:4, 0] = vessel_initial.pos
    params.x_initial[4:7, 0] = vessel_initial.vel
    params.x_initial[7:11, 0] = vessel_initial.rotation
    params.x_initial[11:14, 0] = vessel_initial.omega
        
    #final
    params.x_final[0:1, 0] = vessel_profile.m_dry
    params.x_final[1:4, 0] = vessel_final.pos
    params.x_final[4:7, 0] = vessel_final.vel
    params.x_final[7:11, 0] = vessel_final.rotation
    params.x_final[11:14, 0] = vessel_final.omega
    
    # time of flight guess
    params.s_last = vessel_profile.time_guess
    
    # solver options 凸优化目标函数里各部分权重
    #params.w_delta = solver_options.w_delta
    params.w_nu = solver_options.w_nu
    params.w_delta_s = solver_options.w_delta_s

    integrator = Integrator(vessel_profile, params_super, use_c=False) #use_c is WIP don't use
    
    #  initial guess with GFOLD
    if (verbose):
        print('solving GFOLD for initial guess')
    x_guess, u_guess, m_guess = None, None, None
    for i in range(10): # gfold until feasible
        x_guess, u_guess, m_guess = GFOLD_solver.solve(vessel_profile, vessel_initial, vessel_final, use_c=use_c, verbose=verbose)
        if type(x_guess) == type(None):
            vessel_profile.time_guess *= 1.2
            if (verbose):
                print('GFOLD retry %s, tf=%s' % (i+1, vessel_profile.time_guess))
        else:
            break
    k_space = np.linspace(0, K-1, GFOLD_params.SuperParams().N) # for interpolation in case N != K
    params.s_last = vessel_profile.time_guess
    for k in range(K):

#        alpha1 = (K - k) / K
#        alpha2 = (k / K)
        #xk = params.x_initial * alpha1 + params.x_final * alpha2
        
        xk = np.zeros((14, 1))
        xk[7:11, 0] = np.array([1.0, 0.0, 0.0, 0.0])
        xk[0, 0] = np.interp(k, k_space, m_guess[0, :])#m_guess[0, k]
        for l_i in range(1, 7):
            xk[l_i, 0] = np.interp(k, k_space, x_guess[l_i-1, :])
        uk = np.zeros((3, 1))
        for l_i in range(3):
            uk[l_i, 0] = np.linalg.norm(np.interp(k, k_space, u_guess[l_i, :]))
        params.x_last[:, k] = xk[:, 0]
        #params.u_last[:, k] = xk[0, 0] * -vessel_profile.g  # hover
        #params.u_last_dir[:, k] = -vessel_profile.g / np.linalg.norm(vessel_profile.g)  # thrust dir down
        params.u_last[:, k] = np.array([np.linalg.norm(uk), 0, 0])
        params.u_last_dir[:, k] = params.u_last[:, k] / np.linalg.norm(params.u_last[:, k])
    x_res, u_res, s_res = params.x_last, params.u_last, params.s_last
    
    #迭代求解
    for iteration in range(iterations):
        if (verbose):
            print("Iteration", iteration + 1)

        #（1）积分  准备每个时间步的矩阵A B C S z  过程中需要计算积分，耗时较长
        start_time = time()
        for k in range(0, K - 1):
            params.A[k][:, :], params.B[k][:, :], params.C[k][:, :], params.S[k][:, 0], params.z[k][:, 0] = \
                integrator.solve(params.x_last[:, k], params.u_last[:, k], params.u_last[:, k+1], params.s_last)
        #积分耗时
        if (verbose): 
            print(time() - start_time, "sec to integrate")
        

        # 令w_delta随着迭代数改变
        if isinstance(solver_options.w_delta, type(lambda:0)):
            params.w_delta = solver_options.w_delta(iteration)
        else:
            params.w_delta = solver_options.w_delta
         
#        if solver_options.force_converge:
#            if iteration < solver_options.force_converge_start:
#                params.w_delta = solver_options.w_delta
#            else:
#                params.w_delta = solver_options.w_delta * solver_options.force_converge_amount # w_delta变大以促进delta收敛
#        else:
#            params.w_delta = solver_options.w_delta
        
        #（2）求解凸优化子问题
        start_time = time()
        res, x_res, u_res, s_res, nu_res, delta_res, delta_s_res = solver.solve(params, params_super)
        #凸优化耗时
        if (verbose): 
            print(time() - start_time, "sec to solve")
        
        # 这次迭代的结果作为下一次迭代的输入
        params.x_last = x_res
        params.u_last = u_res
        params.s_last = s_res
        for k in range(K): #计算u的方向的单位向量
            params.u_last_dir[:, k] = params.u_last[:, k] / np.linalg.norm(params.u_last[:, k])
        
        # 收敛情况的评价指标
        delta_norm = np.linalg.norm(delta_res)
        nu_norm = np.linalg.norm(nu_res, ord=1)
        if (verbose):
            print("Flight time:", s_res, end = ' | ')
            print("Delta_norm:", delta_norm, end = ' | ')
            print("Nu_norm:", nu_norm)
        
        #终止条件
        if delta_norm < solver_options.delta_tol and nu_norm < solver_options.nu_tol:
            if (verbose):
                print("Converged after", iteration + 1, "iterations!")
            break
    
    
    #inv normalization
    x_res[0,   :] *= Um
    x_res[1:4, :] *= Ul
    x_res[4:7, :] *= Ul / Ut
    x_res[11:14, 0] *= 1. / Ut
    u_res *= Um * Ul / Ut**2
    s_res *= Ut
    #返回 状态，控制，tf
    return x_res, u_res, s_res


if __name__ == '__main__':
    use_c = False
    if (len(sys.argv) > 1):
        if ('codegen' in sys.argv[1]):
            use_c = True
    
    print('以示例数据进行测试')
    print('模式：%s' % ('cvxpy_codegen生成代码求解' if use_c else 'cvxpy直接求解'))
    x, u, tf = solve(SC_params.VesselProfile.get_default(), 
        SC_params.VesselState.get_default_initial(), 
        SC_params.VesselState.get_default_final(),
        use_c=use_c, verbose=True)
    
    pickle.dump(x.T, open("trajectory/X.p", "wb"))
    pickle.dump(u.T, open("trajectory/U.p", "wb"))

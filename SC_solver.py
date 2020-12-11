

from scipy.integrate import odeint
import numpy as np
import math
from time import time
import pickle, sys

import SC_params
from dynamics_functions_buffer import Dynamics


def solve(vessel_profile, vessel_initial, vessel_final, solver_options=None, params_super=None, use_c = False, verbose = False):
    #default
    if (params_super == None):
        params_super = SC_params.SuperParams()
    if (solver_options == None):
        solver_options = SC_params.SolverOptions()
    
    #solver
    if (use_c):
        import SC_subproblem_gen as solver
    else:
        import SC_subproblem as solver
    
    #successive
    K = params_super.K
    dt = 1 / (K - 1)
    iterations = params_super.iterations
    if (verbose):
        print('K=%s, iterations=%s' % (K, iterations))
    
    params = SC_params.Params(K)
    #sparse
    params.m_dry = vessel_profile.m_dry
    params.cot_gamma_gs = 1 / math.tan(vessel_profile.gamma_gs)
    params.cos_theta_max = math.cos(vessel_profile.theta_max)
    params.omega_max = vessel_profile.omega_max
    params.sec_delta_max = 1 / math.cos(vessel_profile.delta_max)
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
    
    # get matrix functions
    matrix_functions = Dynamics()
    matrix_functions.set_parameters({
        'alpha': 1.0 / 9.806 / vessel_profile.isp,
        'g_I'  : vessel_profile.g,
        'rTB'  : vessel_profile.r_T_B,
        'J'    : vessel_profile.J_B_I,
    })
    
    #  initial guess
    for k in range(K):

        alpha1 = (K - k) / K
        alpha2 = (k / K)
        
        #x_k 和 x_k+1 两个时间点之间的线性插值
        xk = params.x_initial * alpha1 + params.x_final * alpha2
        xk[7:11, 0] = np.array([1.0, 0.0, 0.0, 0.0])
        #print(xk)
        params.x_last[:, k] = xk[:, 0]
        params.u_last[:, k] = xk[0, 0] * -vessel_profile.g  # hover
        params.u_last_dir[:, k] = -vessel_profile.g / np.linalg.norm(vessel_profile.g)  # thrust dir down
    
    #迭代求解
    for iteration in range(iterations):
        if (verbose):
            print("Iteration", iteration + 1)

        start_time = time()
        
        #准备每个时间步的矩阵A B C S z  过程中需要计算积分，耗时较长
        for k in range(0, K - 1):

            # find A_bar, B_bar, C_bar, Sigma_bar = S_bar, and z_bar

            V0 = np.zeros( (idx[-1],) )
            V0[0      : idx[0]] = params.x_last[:, k]              # X at initial step
            V0[idx[0] : idx[1]] = np.eye(14).reshape(-1)  # PhiA at initial step
            
            #积分 V为xABCSZ拼合向量
            V = odeint(
                    ode_dVdt,  # dV/dt function
                    V0,        # initial value of the V container vector
                    (0, dt),   # integrate over the time width of a step
                    args = (params.u_last[:, k], params.u_last[:, k+1], params.s_last, dt,
                        (matrix_functions.A, matrix_functions.B, matrix_functions.f) )
            )

            V = np.array(V)[1, :]
            
            #print(V)

            params.A[k][:, :] = V[idx[0] : idx[1]].reshape((14, 14))
            params.B[k][:, :] = V[idx[1] : idx[2]].reshape((14, 3))
            params.C[k][:, :] = V[idx[2] : idx[3]].reshape((14, 3))
            params.S[k][:, 0] = V[idx[3] : idx[4]]
            params.z[k][:, 0] = V[idx[4] : idx[5]]
        
        if (verbose): #积分耗时
            print(time() - start_time, "sec to integrate")

        #改变w_delta
        if solver_options.force_converge:
            if iteration < solver_options.force_converge_start:
                params.w_delta = solver_options.w_delta
            else:
                params.w_delta = solver_options.w_delta * solver_options.force_converge_amount # w_delta变大以促进delta收敛
        else:
            params.w_delta = solver_options.w_delta

        params.w_nu = solver_options.w_nu
        params.w_delta_s = solver_options.w_delta_s
        
        
        start_time = time()
        
        #solve
        res, x_res, u_res, s_res, nu_res, delta_res, delta_s_res = solver.solve(params, params_super)
        
        if (verbose): #凸优化耗时
            print(time() - start_time, "sec to solve")
        
        # 这次迭代的结果作为下一次迭代的输入
        params.x_last = x_res
        params.u_last = u_res
        params.s_last = s_res
        for k in range(K):
            #print((params.u_last[:, k]).shape, np.linalg.norm(params.u_last[:, k]))
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
    
    # 状态，控制，tf
    return x_res, u_res, s_res


""" The start/end indexes of each variable in the vector V = XABCSz """
idx  = [ 14 ]                # end of x (14,1)
idx += [idx[0] + (14 * 14)]  # end of A (14,14)
idx += [idx[1] + (14 * 3)]   # end of B (14,3)
idx += [idx[2] + (14 * 3)]   # end of C (14,3)
idx += [idx[3] + (14 * 1)]   # end of S (14,1)
idx += [idx[4] + (14 * 1)]   # end of z (14,1)

dV_dt = np.zeros( (idx[-1],) ) #buffer
def ode_dVdt(V, t, u_t, u_t1, sigma, dt, mat_funcs):
    ''' integrate the problem vector, which is defined as:

        V = [x(14),
                Phi_A(14x14), B_bar(14x3), C_bar(14x3), S_bar(14), z_bar(14)]

        V has no correlation to v (velocity) except that it contains v inside
    '''
    #print('ode_dVdt')
    A, B, f = mat_funcs
    
    #dV_dt = np.zeros( (idx[-1],) )
    dV_dt.fill(0)
    
    alpha = (t / dt)
    beta = (1 - t / dt)

    # Gather x and u
    x = V[0 : idx[0]]
    u = u_t + alpha * (u_t1 - u_t)

    ''' Since PhiA(t) * A(t) = Phi(t_end), PhiA(t) = Phi(t_end) * A.inverse '''
    Phi_A_xi = np.linalg.inv(V[idx[0] : idx[1]].reshape((14, 14)))
    Phi      = V[idx[0] : idx[1]].reshape((14, 14))

    z_t  = - np.matmul(A(x, u, sigma), x)
    z_t += - np.matmul(B(x, u, sigma), u)

    dV_dt[0      : idx[0]] = sigma * f(x, u)
    dV_dt[idx[0] : idx[1]] = np.matmul(A(x, u, sigma), Phi).reshape(-1)
    dV_dt[idx[1] : idx[2]] = np.matmul(Phi_A_xi, B(x, u, sigma)).reshape(-1) * alpha
    dV_dt[idx[2] : idx[3]] = np.matmul(Phi_A_xi, B(x, u, sigma)).reshape(-1) * beta
    dV_dt[idx[3] : idx[4]] = np.matmul(Phi_A_xi, f(x, u))
    dV_dt[idx[4] : idx[5]] = np.matmul(Phi_A_xi, z_t)
    
    #print(dV_dt)
    return dV_dt


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

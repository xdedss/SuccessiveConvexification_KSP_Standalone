

from scipy.integrate import odeint
import numpy as np
from old.dynamics_functions import Dynamics
#import sympy

import SC_params


idx  = [ 14 ]                # end of x (14,1)
idx += [idx[0] + (14 * 14)]  # end of A (14,14)
idx += [idx[1] + (14 * 3)]   # end of B (14,3)
idx += [idx[2] + (14 * 3)]   # end of C (14,3)
idx += [idx[3] + (14 * 1)]   # end of S (14,1)
idx += [idx[4] + (14 * 1)]   # end of z (14,1)

class Integrator:
    
    def __init__(self, vessel_profile, params_super):
        self.matrix_functions = Dynamics()
        self.matrix_functions.set_parameters({
            'alpha': 1.0 / 9.806 / vessel_profile.isp,
            'g_I'  : vessel_profile.g,
            'rTB'  : vessel_profile.r_T_B,
            'J'    : vessel_profile.J_B_I,
        })
        self.dt = 1 / (params_super.K - 1)
        
#        x = sympy.MatrixSymbol('x', 14, 1)
#        Phi_A = sympy.MatrixSymbol('Phi_A', 14, 14)
#        B_bar = sympy.MatrixSymbol('B_bar', 14, 3)
#        C_bar = sympy.MatrixSymbol('C_bar', 14, 3)
#        S_bar = sympy.MatrixSymbol('S_bar', 14, 1)
#        z_bar = sympy.MatrixSymbol('z_bar', 14, 1)
#        
#        t = sympy.Symbol('t')
#        u_t = sympy.MatrixSymbol('u_t', 3, 1)
#        u_t1 = sympy.MatrixSymbol('u_t1', 3, 1)
#        s = sympy.Symbol('s')
#        
#        self.ode_dVdt_expr = ode_dVdt_sep(x, Phi_A, B_bar, C_bar, S_bar, z_bar, t, u_t, u_t1, s, self.dt, (self.matrix_functions.A, self.matrix_functions.B, self.matrix_functions.f))

    def solve(self, x_last_k, u_last_k, u_last_k1, s_last):

        V0 = np.zeros( (idx[-1],) )
        V0[0      : idx[0]] = x_last_k            # X at initial step
        V0[idx[0] : idx[1]] = np.eye(14).reshape(-1)  # PhiA at initial step
        
        #积分 V为xABCSZ拼合向量
        V = odeint(
                ode_dVdt,  # dV/dt function
                V0,        # initial value of the V container vector
                (0, self.dt),   # integrate over the time width of a step
                args = (u_last_k, u_last_k1, s_last, self.dt,
                    (self.matrix_functions.A, self.matrix_functions.B, self.matrix_functions.f) )
        )

        V = np.array(V)[1, :]
        
        #print(V)
        
        return (V[idx[0] : idx[1]].reshape((14, 14)),
            V[idx[1] : idx[2]].reshape((14, 3)),
            V[idx[2] : idx[3]].reshape((14, 3)),
            V[idx[3] : idx[4]],
            V[idx[4] : idx[5]] )
                        




def ode_dVdt(V, t, u_t, u_t1, sigma, dt, mat_funcs):
    ''' integrate the problem vector, which is defined as:

        V = [x(14),
                Phi_A(14x14), B_bar(14x3), C_bar(14x3), S_bar(14), z_bar(14)]

        V has no correlation to v (velocity) except that it contains v inside
    '''
    #print('ode_dVdt')
    A, B, f = mat_funcs
    
    dV_dt = np.zeros( (idx[-1],) )
    #dV_dt.fill(0)
    
    alpha = (t / dt)
    beta = (1 - t / dt)

    # Gather x and u
    x = V[0 : idx[0]]
    u = u_t + alpha * (u_t1 - u_t)

    ''' Since PhiA(t) * A(t) = Phi(t_end), PhiA(t) = Phi(t_end) * A.inverse '''
    Phi      = V[idx[0] : idx[1]].reshape((14, 14))
    Phi_A_xi = np.linalg.inv(Phi)
    
    A_xus = A(x, u, sigma)
    B_xus = B(x, u, sigma)
    f_xu = f(x, u)
    
    z_t  = - np.matmul(A_xus, x)
    z_t += - np.matmul(B_xus, u)

    mul_flat = np.matmul(Phi_A_xi, B_xus).reshape(-1)
    dV_dt[0      : idx[0]] = sigma * f_xu
    dV_dt[idx[0] : idx[1]] = np.matmul(A_xus, Phi).reshape(-1)
    dV_dt[idx[1] : idx[2]] = mul_flat * alpha
    dV_dt[idx[2] : idx[3]] = mul_flat * beta
    dV_dt[idx[3] : idx[4]] = np.matmul(Phi_A_xi, f_xu)
    dV_dt[idx[4] : idx[5]] = np.matmul(Phi_A_xi, z_t)
    
    #print(dV_dt)
    return dV_dt
#
#def ode_dVdt_sep(x, Phi_A, B_bar, C_bar, S_bar, z_bar, t, u_t, u_t1, sigma, dt, mat_funcs):
#    A, B, f = mat_funcs
#    dV_dt = sympy.zeros( idx[-1], 1 )
#    
#    alpha = (t / dt)
#    beta = (1 - t / dt)
#    u = u_t + alpha * (u_t1 - u_t)
#    
#    Phi_A_xi = Phi_A.inv()
#    
#    A_xus = A(x, u, sigma)
#    B_xus = B(x, u, sigma)
#    f_xu = f(x, u)
#    
#    z_t  = - (A_xus * x)
#    z_t += - (B_xus * u)
#    
#    mul_PhiB = Phi_A_xi* B_xus
#    return (sigma * f_xu, 
#        A_xus * Phi_A, 
#        mul_PhiB * alpha, 
#        mul_PhiB * beta, 
#        Phi_A_xi * f_xu, 
#        Phi_A_xi * z_t)
#
## sympy 
#def ode_dVdt_s(V, t, u_t, u_t1, sigma, dt, mat_funcs):
#    ''' integrate the problem vector, which is defined as:
#
#        V = [x(14),
#                Phi_A(14x14), B_bar(14x3), C_bar(14x3), S_bar(14), z_bar(14)]
#
#        V has no correlation to v (velocity) except that it contains v inside
#    '''
#    #print('ode_dVdt')
#    A, B, f = mat_funcs
#    
#    dV_dt = sympy.zeros( idx[-1], 1 )
#    #dV_dt.fill(0)
#    
#    alpha = (t / dt)
#    beta = (1 - t / dt)
#    #print(V)
#    # Gather x and u
#    x = V[0 : idx[0]]
#    u = u_t + alpha * (u_t1 - u_t)
#
#    ''' Since PhiA(t) * A(t) = Phi(t_end), PhiA(t) = Phi(t_end) * A.inverse '''
#    Phi      = sympy.Matrix(V[idx[0] : idx[1]]).reshape(14, 14)
#    #Phi_A_xi = Phi.inv()
#    Phi_A_xi = sympy.MatrixSymbol('Phi_A_xi', 14, 14)
#    
#    A_xus = A(x, u, sigma)
#    B_xus = B(x, u, sigma)
#    f_xu = f(x, u)
#    
#    z_t  = - (A_xus* x)
#    z_t += - (B_xus* u)
#    
#    #print((sigma * f_xu).shape)
#    #mul_flat = sympy.Matrix(Phi_A_xi* B_xus).reshape(14 * 3, 1)
#    dV_dt[0      : idx[0], :] = sigma * f_xu
#    dV_dt[idx[0] : idx[1], :] = sympy.Matrix(A_xus* Phi).reshape(14 * 14, 1)
#    dV_dt[idx[1] : idx[2], :] = mul_flat * alpha
#    dV_dt[idx[2] : idx[3], :] = mul_flat * beta
#    dV_dt[idx[3] : idx[4], :] = (Phi_A_xi* f_xu)
#    dV_dt[idx[4] : idx[5], :] = (Phi_A_xi* z_t)
#    
#    #print(dV_dt)
#    return dV_dt



if __name__ == '__main__':
    #testing
    #from dynamics_functions_sympy import Dynamics
    integrator = Integrator(SC_params.VesselProfile.get_default(), SC_params.SuperParams())
    
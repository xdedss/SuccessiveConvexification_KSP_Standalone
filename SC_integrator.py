

from scipy.integrate import odeint
import numpy as np
from old.dynamics_functions import Dynamics
import sympy
import sys

import SC_params


idx  = [ 14 ]                # end of x (14,1)
idx += [idx[0] + (14 * 14)]  # end of A (14,14)
idx += [idx[1] + (14 * 3)]   # end of B (14,3)
idx += [idx[2] + (14 * 3)]   # end of C (14,3)
idx += [idx[3] + (14 * 1)]   # end of S (14,1)
idx += [idx[4] + (14 * 1)]   # end of z (14,1)

class Integrator:
    
    def __init__(self, vessel_profile, params_super, use_c=False):
        self.use_c = use_c
        if (use_c):
            self.parameters = {
                'alpha': 1.0 / 9.806 / vessel_profile.isp,
                'g_I'  : vessel_profile.g,
                'rTB'  : vessel_profile.r_T_B,
                'J'    : vessel_profile.J_B_I,
                'airfric_k' : vessel_profile.airfric_k,
            }
        
        self.matrix_functions = Dynamics()
        self.matrix_functions.set_parameters({
            'alpha': 1.0 / 9.806 / vessel_profile.isp,
            'g_I'  : vessel_profile.g,
            'rTB'  : vessel_profile.r_T_B,
            'J'    : vessel_profile.J_B_I,
            'airfric_k' : vessel_profile.airfric_k,
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
        if (self.use_c):
            V = odeint(
                    ode_dVdt_use_c,  # dV/dt function
                    V0,        # initial value of the V container vector
                    (0, self.dt),   # integrate over the time width of a step
                    args = (u_last_k, u_last_k1, s_last, self.dt,
                        (self.matrix_functions.A, self.matrix_functions.B, self.matrix_functions.f), self.parameters )
            )
        else:
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

def evalMatExpr(expr, **value_dict):
    sym_dict = {}
    for sym in expr.free_symbols:
        sym_dict[sym] = sympy.Matrix(value_dict[sym.name]) if type(value_dict[sym.name])==np.ndarray else value_dict[sym.name]
    s = expr.subs(sym_dict)
    #print(s.free_symbols)
    return np.array(s.doit()).astype('float')
#np.array((m*2).subs({m:sympy.Matrix(np.array([[1,2],[3,4]]))}).doit(), dtype=float)

# buggy dont use
def ode_dVdt_use_c(V, t, u_t, u_t1, sigma, dt, mat_funcs, parameters):
#    from sc_integrator_codegen.wrapper_module_0 import autofunc_c as func_res_x
#    from sc_integrator_codegen.wrapper_module_1 import autofunc_c as func_res_phi
#    from sc_integrator_codegen.wrapper_module_2 import autofunc_c as func_res_b
#    from sc_integrator_codegen.wrapper_module_3 import autofunc_c as func_res_c
#    from sc_integrator_codegen.wrapper_module_4 import autofunc_c as func_res_s
#    from sc_integrator_codegen.wrapper_module_5 import autofunc_c as func_res_z
    alpha = parameters['alpha']
    g_I = parameters['g_I'].reshape((3, 1))
    rTB = parameters['rTB'].reshape((3, 1))
    J = parameters['J']
    
    V = V.reshape((idx[-1], 1))
    u_t = u_t.reshape((3, 1))
    u_t1 = u_t1.reshape((3, 1))

    ''' Since PhiA(t) * A(t) = Phi(t_end), PhiA(t) = Phi(t_end) * A.inverse '''
    Phi      = V[idx[0] : idx[1]].reshape((14, 14))
    Phi_A_xi = np.linalg.inv(Phi)
    
    dV_dt = np.zeros( (idx[-1],) )
    
#    res_x = func_res_x(J, V, alpha, dt, g_I, rTB, sigma, t, u_t, u_t1)
#    res_phi = func_res_phi(J, Phi, V, dt, sigma, t, u_t, u_t1)
#    res_b = func_res_b(J, Phi_A_xi, V, alpha, dt, rTB, sigma, t, u_t, u_t1)
#    res_c = func_res_c(J, Phi_A_xi, V, alpha, dt, rTB, sigma, t, u_t, u_t1)
#    res_s = func_res_s(J, Phi_A_xi, V, alpha, dt, g_I, rTB, t, u_t, u_t1)
#    res_z = func_res_z(J, Phi_A_xi, V, alpha, dt, rTB, sigma, t, u_t, u_t1)
    
#    res_x = func_res_x(J=J, V=V, alpha=alpha, dt=dt, g_I=g_I, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_phi = func_res_phi(J=J, Phi=Phi, V=V, dt=dt, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_b = func_res_b(J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_c = func_res_c(J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_s = func_res_s(J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, g_I=g_I, rTB=rTB, t=t, u_t=u_t, u_t1=u_t1)
#    res_z = func_res_z(J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
    
    # test: direct eval
    res_x = evalMatExpr(func_res_x, J=J, V=V, alpha=alpha, dt=dt, g_I=g_I, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_phi = evalMatExpr(func_res_phi, J=J, Phi=Phi, V=V, dt=dt, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_b = evalMatExpr(func_res_b, J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_c = evalMatExpr(func_res_c, J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
#    res_s = evalMatExpr(func_res_s, J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, g_I=g_I, rTB=rTB, t=t, u_t=u_t, u_t1=u_t1)
#    res_z = evalMatExpr(func_res_z, J=J, Phi_A_xi=Phi_A_xi, V=V, alpha=alpha, dt=dt, rTB=rTB, sigma=sigma, t=t, u_t=u_t, u_t1=u_t1)
    
    # test: original
    A, B, f = mat_funcs
    alpha = (t / dt)
    beta = (1 - t / dt)
    x = V[0 : idx[0]]
    u = u_t + alpha * (u_t1 - u_t)
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
    dV_dt[idx[3] : idx[4]] = np.matmul(Phi_A_xi, f_xu).reshape(-1)
    dV_dt[idx[4] : idx[5]] = np.matmul(Phi_A_xi, z_t).reshape(-1)
    
#    dV_dt[0      : idx[0]] = res_x.reshape(-1)
#    dV_dt[idx[0] : idx[1]] = res_phi.reshape(-1)
#    dV_dt[idx[1] : idx[2]] = res_b.reshape(-1)
#    dV_dt[idx[2] : idx[3]] = res_c.reshape(-1)
#    dV_dt[idx[3] : idx[4]] = res_s.reshape(-1)
#    dV_dt[idx[4] : idx[5]] = res_z.reshape(-1)
    
    #print(dV_dt)
    return dV_dt

# sympy 
def ode_dVdt_expr():
    ''' integrate the problem vector, which is defined as:

        V = [x(14),
                Phi_A(14x14), B_bar(14x3), C_bar(14x3), S_bar(14), z_bar(14)]

        V has no correlation to v (velocity) except that it contains v inside
    '''
    
    import sympy
    from old.dynamics_functions_sympy import Dynamics as Dynamics_sympy
    
    matrix_functions = Dynamics_sympy()
    matrix_functions.set_parameters({
        'alpha': sympy.Symbol('alpha'),
        'g_I'  : sympy.MatrixSymbol('g_I', 3, 1),
        'rTB'  : sympy.MatrixSymbol('rTB', 3, 1),
        'J'    : sympy.MatrixSymbol('J', 3, 3),
    })
    
    V = sympy.MatrixSymbol('V', idx[-1], 1)
    t = sympy.Symbol('t')
    u_t = sympy.MatrixSymbol('u_t', 3, 1)
    u_t1 = sympy.MatrixSymbol('u_t1', 3, 1)
    sigma = sympy.Symbol('sigma')
    dt = sympy.Symbol('dt')
    
    A, B, f = (matrix_functions.A, matrix_functions.B, matrix_functions.f)
    
    alpha = (t / dt)
    beta = (1 - t / dt)
    # Gather x and u
    x = V[0 : idx[0]]
    u = u_t + alpha * (u_t1 - u_t)

    ''' Since PhiA(t) * A(t) = Phi(t_end), PhiA(t) = Phi(t_end) * A.inverse '''
    #Phi      = sympy.Matrix(V[idx[0] : idx[1]]).reshape(14, 14)
    Phi = sympy.MatrixSymbol('Phi', 14, 14)
    #Phi_A_xi = Phi.inv()
    Phi_A_xi = sympy.MatrixSymbol('Phi_A_xi', 14, 14)
    
    A_xus = A(x, u, sigma)
    B_xus = B(x, u, sigma)
    f_xu = f(x, u)
    
    z_t  = - (A_xus* x)
    z_t += - (B_xus* u)
    
    #print((sigma * f_xu).shape)
    #mul_flat = sympy.Matrix(Phi_A_xi* B_xus).reshape(14 * 3, 1)
    mul_ab = Phi_A_xi* B_xus
    res_x = sigma * f_xu
#    res_phi = sympy.Matrix(A_xus* Phi)
#    res_b = sympy.Matrix(mul_ab * alpha)
#    res_c = sympy.Matrix(mul_ab * beta)
#    res_s= sympy.Matrix(Phi_A_xi* f_xu)
#    res_z = sympy.Matrix(Phi_A_xi* z_t)
    
    res_phi = (A_xus* Phi)
    res_b = (mul_ab * alpha)
    res_c = (mul_ab * beta)
    res_s= (Phi_A_xi* f_xu)
    res_z = (Phi_A_xi* z_t)
    
    #print(dV_dt)
    return (res_x, res_phi, res_b, res_c, res_s, res_z)


#(func_res_x, func_res_phi, func_res_b, func_res_c, func_res_s, func_res_z) = ode_dVdt_expr()
def codegen(codegen_path):
    from sympy.utilities.autowrap import autowrap
    print('codegen...')
    expr = ode_dVdt_expr()
    names = ['res_x', 'res_phi', 'res_b', 'res_c', 'res_s', 'res_z']
    print('generating...')
    for i in range(len(expr)):
        print('free sym of %s:' % names[i])
        print(expr[i].free_symbols)
        autowrap(expr[i], language='C', backend='cython', tempdir=codegen_path)

#    free sym of res_x:
#    {u_t, J, t, sigma, u_t1, alpha, g_I, V, rTB, dt}
#    free sym of res_phi:
#    {u_t, J, t, sigma, u_t1, Phi, V, dt}
#    free sym of res_b:
#    {u_t, J, t, sigma, u_t1, alpha, Phi_A_xi, V, rTB, dt}
#    free sym of res_c:
#    {u_t, J, t, sigma, u_t1, alpha, Phi_A_xi, V, rTB, dt}
#    free sym of res_s:
#    {u_t, J, t, u_t1, alpha, g_I, Phi_A_xi, V, rTB, dt}
#    free sym of res_z:
#    {u_t, J, t, sigma, u_t1, alpha, Phi_A_xi, V, rTB, dt}


if __name__ == '__main__':
    #testing
    #from dynamics_functions_sympy import Dynamics
    integrator = Integrator(SC_params.VesselProfile.get_default(), SC_params.SuperParams())
    
    if (len(sys.argv) > 2 and sys.argv[1] == 'codegen'):
        codegen_path = sys.argv[2]
        codegen(codegen_path)
    else:
        print("invalid input")
        print(sys.argv)


#
#from sympy import *
#from sympy.abc import x, y, z
#from sympy.utilities.autowrap import autowrap
#expr = ((x - y + z)**(13)).expand()
#binary_func = autowrap(expr, language='C', backend='cython', tempdir='sc_integrator_codegen')
#
#
#try:
#    from sc_integrator_codegen.wrapper_module_0 import autofunc_c as odeint_c
#except ImportError:


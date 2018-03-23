'''
    Ackimese & Szmuk : Successive Convexification for 6-DoF Mars Landing

        -> https://arxiv.org/abs/1802.03827

    FORKED FROM Sven Niederberger's implementation
        -> https://github.com/EmbersArc/SuccessiveConvexification

    CHANGES FROM Niederberger's implemenation:
        - Due to intense hatred of variables ending with "_", replaced all
            cvx variables ending with "_" with proper names
            - v instead of _ denotes a cvx Variable
            - parm instead of _ denotes a cvx Parameter

        - Replaced a great number of () with [] to clarify defintions

        - Replaced all dynamics functions in parameters with the sympy
          generated ones in dynamics_functions.py
'''

from scipy.integrate import odeint
from parameters import *
from time import time
import cvxpy as cvx
import pickle


'''
    state variable

               0   1   2   3   4   5   6   7   8   9  10  11  12  13
          x = [m, r0, r1, r2, v0, v1, v2, q0, q1, q2, q3, w0, w1, w2].T

    control variable, thrust

              u = [T0, T1, T2] in body axes
'''

Xi = np.empty(shape=[K, 14])  # axis 0 = discretization step,
Ui = np.empty(shape=[K, 3])   # axis 1 = state / control variable axis

x_init  = np.concatenate([ [m_wet,],
                                    r_I_init, v_I_init, q_B_I_init, w_B_init])
x_final = np.concatenate([ [m_dry,],
                                r_I_final, v_I_final, q_B_I_final, w_B_final])

# SETUP ---------------------------------------------------------------------

# Variables:
Xv       = cvx.Variable((K, 14))
Uv       = cvx.Variable((K, 3) )
nuv      = cvx.Variable((14 * (K - 1)))
delta    = cvx.Variable(K)
sigmav   = cvx.Variable()
delta_sv = cvx.Variable()

# Parameters:
A_bar_parm = cvx.Parameter((K, 14 * 14))
B_bar_parm = cvx.Parameter((K, 14 * 3))
C_bar_parm = cvx.Parameter((K, 14 * 3))
S_bar_parm = cvx.Parameter((K, 14))
z_bar_parm = cvx.Parameter((K, 14))

X_last_parm        = cvx.Parameter((K, 14))
U_last_parm        = cvx.Parameter((K, 3))
s_last_parm        = cvx.Parameter(nonneg=True)
w_delta_parm       = cvx.Parameter(nonneg=True)
w_nu_parm          = cvx.Parameter(nonneg=True)
w_delta_sigma_parm = cvx.Parameter(nonneg=True)

# Boundary conditions:
constraints = [
    Xv[0, 0] == x_init[0],
    Xv[0, 1:4] == x_init[1:4],
    Xv[0, 4:7] == x_init[4:7],
    # Xv[0, 7:11] == x_init[7:11],  # initial attitude is free
    Xv[0, 11:14] == x_init[11:14],

    # Xv[0, 0] == x_final[0], # final mass is free
    Xv[K - 1, 1:4] == x_final[1:4],
    Xv[K - 1, 4:7] == x_final[4:7],
    Xv[K - 1, 7:11] == x_final[7:11],
    Xv[K - 1, 11:14] == x_final[11:14],

    Uv[K - 1, 1] == 0,
    Uv[K - 1, 2] == 0
]

# Dynamics:
for k in range(K - 1):
    constraints += [
        Xv[k + 1, :] ==
        cvx.reshape(A_bar_parm[k, :], (14, 14)) * Xv[k, :]
        + cvx.reshape(B_bar_parm[k, :], (14, 3)) * Uv[k, :]
        + cvx.reshape(C_bar_parm[k, :], (14, 3)) * Uv[k + 1, :]
        + S_bar_parm[k, :] * sigmav
        + z_bar_parm[k, :]
        + nuv[k * 14:(k + 1) * 14]
    ]

# State constraints:
constraints += [Xv[:, 0] >= m_dry]
for k in range(K):
    constraints += [
        tan_gamma_gs * cvx.norm(Xv[k, 2: 4]) <= Xv[k, 1],
        cos_theta_max <= 1 - 2 * cvx.sum_squares(Xv[k, 9:11]),
        cvx.norm(Xv[k, 11: 14]) <= w_B_max
    ]

# Control constraints:
for k in range(K):
    B_g = U_last_parm[k, :] / cvx.norm(U_last_parm[k, :])
    constraints += [
        T_min <= B_g * Uv[k, :],
        cvx.norm(Uv[k, :]) <= T_max,
        cos_delta_max * cvx.norm(Uv[k, :]) <= Uv[k, 0]
    ]

# Trust regions:
for k in range(K):
    dx = Xv[k, :] - X_last_parm[k, :]
    du = Uv[k, :] - U_last_parm[k, :]
    constraints += [
        cvx.sum_squares(dx) + cvx.sum_squares(du) <= delta[k]
    ]

ds = sigmav - s_last_parm
constraints += [cvx.norm(ds, 2) <= delta_sv]

# Objective:
objective = cvx.Minimize(
    sigmav  # minimize flight time
    #-Xv[0,-1]  # minimize fuel expenditure / maximize terminal fuel
        + w_nu_parm     * cvx.norm(nuv, 1)  # virtual control
        + w_delta_parm  * cvx.norm(delta)   # trust region on dynamics
        + w_delta_sigma_parm * cvx.norm(delta_sv, 1)  # trust region on sigma
)

prob = cvx.Problem(objective, constraints)

if not prob.is_dcp():
    raise(Exception("Problem does not follow DCP rules"))

# Initialization ------------------------------------------------------------

sigma = t_f_guess

for k in range(K):

    alpha1 = (K - k) / K
    alpha2 = (k / K)

    # Linear interpolate between start and end values for linearization guess

    mk   = [alpha1 * x_init[0] + alpha2 * x_final[0], ]
    rk =  alpha1 * x_init[1:4] + alpha2 * x_final[1:4]
    vk =  alpha1 * x_init[4:7] + alpha2 * x_final[4:7]
    wk =  alpha1 * x_init[11:14] + alpha2 * x_final[11:14]
    qk = np.array([1.0, 0.0, 0.0, 0.0])

    Xi[k, :] = np.concatenate([mk, rk, vk, qk, wk])
    Ui[k, :] = mk * -g_I  # hover


""" The start/end indexes of each variable in the vector V = XABCSz """
idx  = [ 14 ]                # end of x (14,1)
idx += [idx[0] + (14 * 14)]  # end of A (14,14)
idx += [idx[1] + (14 * 3)]   # end of B (14,3)
idx += [idx[2] + (14 * 3)]   # end of C (14,3)
idx += [idx[3] + (14 * 1)]   # end of S (14,1)
idx += [idx[4] + (14 * 1)]   # end of z (14,1)


def ode_dVdt(V, t, u_t, u_t1, sigma):
    ''' integrate the problem vector, which is defined as:

        V = [x(14),
                Phi_A(14x14), B_bar(14x3), C_bar(14x3), S_bar(14), z_bar(14)]

        V has no correlation to v (velocity) except that it contains v inside
    '''

    dV_dt = np.zeros( (idx[-1],) )

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

    return dV_dt


# Successive Convexification -------------------------------------------------
for iteration in range(iterations):

    print("Iteration", iteration + 1)

    start_time = time()

    A_bar = np.zeros([K, 14, 14])
    B_bar = np.zeros([K, 14, 3])
    C_bar = np.zeros([K, 14, 3])
    S_bar = np.zeros([K, 14])
    z_bar = np.zeros([K, 14])

    for k in range(0, K - 1):

        # find A_bar, B_bar, C_bar, Sigma_bar = S_bar, and z_bar

        V0 = np.zeros( (idx[-1],) )
        V0[0      : idx[0]] = Xi[k, :]                # X at initial step
        V0[idx[0] : idx[1]] = np.eye(14).reshape(-1)  # PhiA at initial step

        V = odeint(
                        ode_dVdt,  # dV/dt function
                        V0,        # initial value of the V container vector
                        (0, dt),   # integrate over the time width of a step
                        args = (Ui[k, :], Ui[k + 1, :], sigma)
                )

        V = np.array(V)[1, :]

        A_bar[k, :, :] = V[idx[0] : idx[1]].reshape((14, 14))
        B_bar[k, :, :] = V[idx[1] : idx[2]].reshape((14, 3))
        C_bar[k, :, :] = V[idx[2] : idx[3]].reshape((14, 3))
        S_bar[k,    :] = V[idx[3] : idx[4]]
        z_bar[k,    :] = V[idx[4] : idx[5]]

    X_last_parm.value = Xi
    U_last_parm.value = Ui
    s_last_parm.value = sigma

    A_bar_parm.value  = A_bar.reshape((K, 14 * 14), order='F')
    B_bar_parm.value  = B_bar.reshape((K, 14 * 3), order='F')
    C_bar_parm.value  = C_bar.reshape((K, 14 * 3), order='F')
    S_bar_parm.value  = S_bar
    z_bar_parm.value  = z_bar

    #  weights are defined and imported from parameters.py

    if force_converge['active']:
        if iteration < force_converge['start']:
            w_delta_parm.value  =  w_delta
        else:
            w_delta_parm.value  =  w_delta * force_converge['amount']
    else:
        w_delta_parm.value = w_delta

    w_nu_parm.value = w_nu
    w_delta_sigma_parm.value = w_delta_sigma

    print("Solving problem. ", time() - start_time, "sec to integrate")

    prob.solve(verbose=True, solver='ECOS')
    # CVX ----------------------------------------------------------------------------------------------------------

    X = Xv.value
    U = Uv.value
    sigma = sigmav.value

    delta_norm = np.linalg.norm(delta.value)
    nu_norm = np.linalg.norm(nuv.value, ord=1)

    print("Flight time:", sigmav.value, end = ' | ')
    print("Delta_norm:", delta_norm, end = ' | ')
    print("Nu_norm:", nu_norm)

    if delta_norm < delta_tol and nu_norm < nu_tol:
        print("Converged after", iteration + 1, "iterations!")
        break

pickle.dump(X, open("trajectory/X.p", "wb"))
pickle.dump(U, open("trajectory/U.p", "wb"))

print('Saved. Program complete!')

from dynamics_functions import Dynamics
import numpy as np

# Trajectory points
K = 70
dt = 1 / (K - 1)

# Max solver iterations
iterations = 15

# increase w_delta for iterations above 'start'
# this puts major pressure on the dynamics trust region to converge
force_converge = {
                    'active'  : False,
                    'start'   : 4,
                    'amount'  : 1e3
                }

# Mass
m_wet = 2.0
m_dry = 1.0

# Flight time guess
t_f_guess = 5.

# # Weight constants
w_nu = 1e5
w_delta = 1e-3
w_delta_sigma = 1e-1

# Exit conditions
nu_tol = 1e-8
delta_tol = 1e-3

# State constraints
r_I_init = np.array([20., 5., 5.])
v_I_init = np.array([-4, -3, -2])
q_B_I_init = np.array([1.0, 0.0, 0.0, 0.0])
w_B_init = np.array([0., 0., 0.])

r_I_final = np.array([0., 0., 0.])
v_I_final = np.array([0., 0., 0.])
q_B_I_final = np.array([1.0, 0.0, 0.0, 0.0])
w_B_final = np.array([0., 0., 0.])

w_B_max = np.deg2rad(60)

# Angles
cos_delta_max = np.cos(np.deg2rad(20))
cos_theta_max = np.cos(np.deg2rad(90))
tan_gamma_gs = np.tan(np.deg2rad(42))

# Angular moment of inertia
J_B_I = 1e-2 * np.eye(3)

# Vector from thrust point to CoM
r_T_B = np.array((-1e-2, 0., 0.))

# Gravity
g_I = np.array((-1., 0., 0.))

# Thrust limits
T_min = 2.0
T_max = 5.0

# Fuel consumption
alpha_m = 0.01

parms = {
            'alpha': alpha_m,
            'g_I'  : g_I,

            'rTB'  : r_T_B,
            'J'    : J_B_I,
}

# Linearized state matrices:
matrix_functions = Dynamics()
matrix_functions.set_parameters(parms)

# Add to local namespace so that the main program can see them here
A = matrix_functions.A
B = matrix_functions.B
f = matrix_functions.f

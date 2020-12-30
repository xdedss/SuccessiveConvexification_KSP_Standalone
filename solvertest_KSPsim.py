
import numpy as np

import SC_solver, SC_params
import trajectory.plot as plot

import pickle


vessel = SC_params.VesselProfile.get_default()
state_initial = SC_params.VesselState.get_default_initial()
state_final = SC_params.VesselState.get_default_final()


#vessel.T_min = 1.0
#vessel.T_max = 5.0
#state_initial.pos = np.array([20., 5., 3.])
#state_initial.vel = np.array([-4., -0., -3.])


vessel.isp = 299.5
vessel.g = np.array((-9.807, 0., 0.)) # gravity
vessel.m_dry = 19770.810546875
vessel.gamma_gs = np.deg2rad(10) # glide slope
vessel.theta_max = np.deg2rad(90) # tilt
vessel.theta_max = np.linspace(np.deg2rad(90), np.deg2rad(10), SC_params.SuperParams().K)
vessel.omega_max = np.deg2rad(60) # rotation vel
vessel.delta_max = np.deg2rad(10) # gimbal
vessel.T_min = 190162.0
vessel.T_max = 760648.0
vessel.r_T_B = np.array([-4.49584946,  0.       ,   0.        ])
vessel.J_B_I = np.array([[ 3.43008688e+05,  9.74064052e-01, -2.31322622e+00],
 [ 9.73825634e-01,  3.87805234e+04,  8.40314293e+00],
 [-2.31325603e+00,  8.40314293e+00,  3.49099750e+05]])
vessel.airfric_k = 14
vessel.time_guess = 25

state_initial.mass = 25583.166015625
state_initial.pos = np.array([1400.        ,  165.28223741 * 0.9,  247.82434027 * 0.9])
state_initial.vel = np.array([-98.40400201 , -6.23042187, -10.52199281])
state_initial.rotation = np.array([ 0.70710678 , 0.         , 0.58785953, -0.39296459])
state_initial.omega = np.array([0., 0., 0.])

state_final.mass = 25583.166015625
state_final.pos = np.array([200.  , 0. , 0.])
state_final.vel = np.array([-89.27567157 ,  0.       ,    0.        ])
state_final.pos = np.array([50.  , 0. , 0.])
state_final.vel = np.array([-45 ,  0.       ,    0.        ])
state_final.rotation = np.array([1. ,0., 0. ,0.])
state_final.omega = np.array([0., 0., 0.])

solver_options = SC_params.SolverOptions()
#solver_options.w_nu = 1e5
#solver_options.w_delta = lambda i:(1e-3 if i < 4 else (1 if i < 8 else 10))
solver_options.w_delta = lambda i:(1e-3 * (2 ** i))
#solver_options.w_delta_s = 1e-1
#solver_options.nu_tol = 1e-8
#solver_options.delta_tol = 1e-3
#solver_options.force_converge = True
#solver_options.force_converge_start = 4
#solver_options.force_converge_amount = 1e3

Ut = vessel.time_guess / 5.
Ul = state_initial.pos[0] / 20.
Um = state_initial.mass / 1.5
vessel = SC_params.normalize(vessel, Ut, Ul, Um)
state_initial = SC_params.normalize(state_initial, Ut, Ul, Um)
state_final = SC_params.normalize(state_final, Ut, Ul, Um)

print(vessel.__dict__)
print(state_initial.__dict__)

x, u, tf = SC_solver.solve(vessel, state_initial, state_final, 
    solver_options=solver_options, use_c=True, verbose=True)


plot.plot_3D(x.T, u.T)
#plot.plot_X(x.T, u.T)

pickle.dump(x.T, open("trajectory/X.p", "wb"))
pickle.dump(u.T, open("trajectory/U.p", "wb"))



# Todo：自动normalize，测试真实尺度数据
# Todo: 用GFOLD初始化

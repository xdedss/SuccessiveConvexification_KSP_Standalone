
import numpy as np

import SC_solver, SC_params
import GFOLD_solver, GFOLD_params
import trajectory.plot as plot

import pickle

N=GFOLD_params.SuperParams().N
K=SC_params.SuperParams().K

vessel = SC_params.VesselProfile.get_default()
state_initial = SC_params.VesselState.get_default_initial()
state_final = SC_params.VesselState.get_default_final()


#vessel.T_min = 1.0
#vessel.T_max = 5.0
#state_initial.pos = np.array([20., 5., 3.])
#state_initial.vel = np.array([-4., -0., -3.])


thrust = 24000
vessel.isp = 203.94
vessel.g = np.array((-3.71, 0., 0.)) # gravity
vessel.m_dry = 2.0e3
vessel.gamma_gs = np.deg2rad(30) # glide slope
vessel.theta_max = np.deg2rad(60) # tilt
vessel.theta_max = np.linspace(vessel.theta_max, np.deg2rad(10), K)
vessel.omega_max = np.deg2rad(60) # rotation vel
vessel.delta_max = np.deg2rad(20) # gimbal
vessel.T_min = 0.2*thrust
vessel.T_max = 0.8*thrust
vessel.r_T_B = np.array((-0.5, 0., 0.))
vessel.J_B_I = 20000 * np.eye(3)
vessel.airfric_k = 5
vessel.time_guess = 50

state_initial.mass = 2.3e3
state_initial.pos = np.array([1400., 450., -330.])
state_initial.vel = np.array([-20., 40., 40.])
state_initial.rotation = np.array([1.0, 0.0, 0.0, 0.0])
state_initial.omega = np.array([0., 0., 0.])

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
#vessel = SC_params.normalize(vessel, Ut, Ul, Um)
#state_initial = SC_params.normalize(state_initial, Ut, Ul, Um)

print(vessel.__dict__)
print(state_initial.__dict__)


#x, u, m = GFOLD_solver.solve(vessel, state_initial, state_final, use_c=False, verbose=True)
#x = np.vstack([m, x, np.ones((1, N)), np.zeros((6, N))])
#u = u * m


x, u, tf = SC_solver.solve(vessel, state_initial, state_final, 
    solver_options=solver_options, use_c=True, verbose=True)

x[0,   :] /= Um
x[1:4, :] /= Ul
x[4:7, :] /= Ul / Ut
x[11:14, 0] /= 1. / Ut
u /= Um * Ul / Ut**2
tf /= Ut



plot.plot_3D(x.T, u.T)
#plot.plot_X(x.T, u.T)

pickle.dump(x.T, open("trajectory/X.p", "wb"))
pickle.dump(u.T, open("trajectory/U.p", "wb"))



# Todo：自动normalize，测试真实尺度数据
# Todo: 用GFOLD初始化

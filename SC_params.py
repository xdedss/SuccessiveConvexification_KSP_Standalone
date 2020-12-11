
import numpy as np

class SuperParams:
    def __init__(self):
        self.K = 70
        self.iterations = 15

class SolverOptions:
    def __init__(self):
        self.w_nu = 1e5
        self.w_delta = 1e-3
        self.w_delta_s = 1e-1
        #exit condition
        self.nu_tol = 1e-8
        self.delta_tol = 1e-3

class Params:
    def __init__(self, K): #K为时间步数
        self.A = np.zeros((K, 14, 14))
        self.B = np.zeros((K, 14, 3))
        self.C = np.zeros((K, 14, 3))
        self.S = np.zeros((K, 14, 1))
        self.z = np.zeros((K, 14, 1))
        self.x_last = np.zeros((14, K))
        self.u_last = np.zeros((3, K))
        self.u_last_dir = np.zeros((3, K))
        self.s_last = None
        self.w_nu = None
        self.w_delta = None
        self.w_delta_s = None
        
        self.x_initial = np.zeros((14, 1))
        self.x_final = np.zeros((14, 1))
        #sparse
        self.m_dry = None
        self.cot_gamma_gs = None
        self.cos_theta_max = None
        self.omega_max = None
        self.sec_delta_max = None
        self.T_max = None
        self.T_min = None
    
    def check(self):
        #print(self.A)
        pass


class VesselProfile:
    def __init__(self):
        self.isp = None
        self.g = np.array([-9.806, 0.0, 0.0]) # gravity
        self.m_dry = None
        self.gamma_gs = None # glide slope
        self.theta_max = None # tilt
        self.omega_max = None # rotation vel
        self.delta_max = None # gimbal
        self.T_min = None
        self.T_max = None
        self.r_T_B = None
        self.J_B_I = None

        self.time_guess = None
    
    def get_default():
        default_vessel = VesselProfile()
        default_vessel.isp = 10.194
        default_vessel.g = np.array((-1., 0., 0.)) # gravity
        default_vessel.m_dry = 1.0
        default_vessel.gamma_gs = np.deg2rad(20) # glide slope
        default_vessel.theta_max = np.deg2rad(90) # tilt
        default_vessel.omega_max = np.deg2rad(60) # rotation vel
        default_vessel.delta_max = np.deg2rad(20) # gimbal
        default_vessel.T_min = 2.0
        default_vessel.T_max = 5.0
        default_vessel.r_T_B = np.array((-1e-2, 0., 0.))
        default_vessel.J_B_I = 1e-2 * np.eye(3)
        default_vessel.time_guess = 5.
        return default_vessel
        

class VesselState:
    def __init__(self):
        self.mass = None
        self.pos = None
        self.vel = None
        self.rotation = None
        self.omega = None
    
    def get_default_initial():
        default_initial = VesselState()
        default_initial.mass = 2.0
        default_initial.pos = np.array([20., 5., 5.])
        default_initial.vel = np.array([-4, -3, -2])
        default_initial.rotation = np.array([1.0, 0.0, 0.0, 0.0])
        default_initial.omega = np.array([0., 0., 0.])
        return default_initial
    
    def get_default_final():
        default_final = VesselState()
        default_final.pos = np.array([0., 0., 0.])
        default_final.vel = np.array([0., 0., 0.])
        default_final.rotation = np.array([1.0, 0.0, 0.0, 0.0])
        default_final.omega = np.array([0., 0., 0.])
        return default_final


class SuperParams:
    def __init__(self):
        self.K = 10


class Params:
    def __init__(self):
        self.A = None
        self.B = None
        self.C = None
        self.S = None
        self.z = None
        self.x_last = None
        self.u_last = None
        self.u_last_dir = None
        self.s_last = None
        self.w_nu = None
        self.w_delta = None
        self.w_delta_s = None
        
        self.x_initial = None
        self.x_final = None
        #sparse
        self.m_dry = None
        self.cot_gamma_gs = None
        self.cos_theta_max = None
        self.omega_max = None
        self.sec_delta_max = None
        self.T_max = None
        self.T_min = None




# 空间运算用

import math
import numpy as np
import numpy.linalg as npl

def clamp(num, maxnum, minnum):
    if maxnum < minnum:
        return clamp(num, minnum, maxnum)
    if num > maxnum:
        return maxnum
    elif num < minnum:
        return minnum
    return num

def clamp_mag(vec, maxmag):
    mag = npl.norm(vec)
    if mag > maxmag:
        return vec / mag * maxmag
    return vec

def lerp(vec1, vec2, t):
    return t * vec2 + (1-t) * vec1

def sgn(f):
    if f > 0:
        return 1
    elif f < 0:
        return -1
    return 0

def v2(f1, f2):
    return np.array((f1, f2), dtype='float')

def v3(f1, f2, f3):
    return np.array((f1, f2, f3), dtype='float')

def vec(*args):
    if (len(args) == 1):
        args = args[0]
    return np.array(args, dtype='float')

def normalize(vec):
    return vec / npl.norm(vec)

# x y z w 
def quat(axis, angle):
    (x, y, z) = axis
    s = math.sin(angle / 2)
    c = math.cos(angle / 2)
    axis = v3(x+.0, y+.0, z+.0)
    axis /= npl.norm(axis)
    return (s * axis[0], s * axis[1], s * axis[2], c)

# x y z w 
def rotation_mat(q):
    q = q / np.linalg.norm(q)
    (x, y, z, w) = q
    return np.mat([
        [1-2*y**2-2*z**2, 2*x*y+2*w*z, 2*x*z-2*w*y],
        [2*x*y-2*w*z, 1-2*x**2-2*z**2, 2*y*z+2*w*x],
        [2*x*z+2*w*y, 2*y*z-2*w*x, 1-2*x**2-2*y**2]
    ]).T

# 转列向量右乘
def transform(vec1, mat):
    res = (mat * np.mat(vec1).T).T
    return vec(res[0,0], res[0,1], res[0,2])

def angle_around_axis(v1, v2, axis):
    axis = normalize(axis)
    v1 = normalize(np.cross(v1, axis))
    v2 = normalize(np.cross(v2, axis))
    direction = sgn(np.dot(np.cross(v1, v2), axis))
    return direction * math.acos(np.dot(v1, v2))

def norm_deg(angle, center=0):
    angle = angle % 360
    while (angle < center - 180):
        angle += 360
    while (angle >= center + 180):
        angle -= 360
    return angle

class FreeInertialControl:
    # 假设：几乎无干扰
    def __init__(self):
        self.unified_authority = 1
        self.redundancy = 0.1
        self.kp = 1
        self.kd = 0
    
    def update(self, error, velocity):
        acc = self.unified_authority * (1 - self.redundancy)
        target_vel = math.sqrt(2 * max(abs(error) - 0.5 * deg2rad, 0) * acc) * sgn(-error)
        self.result = (target_vel - velocity) * self.kp - velocity * self.kd
        return self.result

class PID:
    def __init__(self):
        self.ep = True
        self.ei = True
        self.ed = True
        self.kp = 1
        self.ki = 0
        self.kd = 1
        self.sd = 0
        self.diff = 0
        self.integral = 0
        self.integral_limit = 1
        self.error_prev = 0
        self.first = True
        self.second = True
        self.dumpf = None
    
    def update(self, error, dt):
        if self.first:
            self.first = False
            self.error_prev = error
        elif self.second:
            self.second = False
            self.diff = (error - self.error_prev) / dt
        
        self.integral += error * dt * self.ki
        self.integral = clamp(self.integral, self.integral_limit, -self.integral_limit)
        self.diff = lerp(self.diff, (error - self.error_prev) / dt, 1-self.sd)
        p = -error * self.kp
        i = -self.integral
        d = -self.diff * self.kd
        self.result = p * (1 if self.ep else 0) + i * (1 if self.ei else 0) + d * (1 if self.ed else 0)
        
        self.p_prev = p
        self.i_prev = i
        self.d_prev = d
        self.error_prev = error
        return self.result




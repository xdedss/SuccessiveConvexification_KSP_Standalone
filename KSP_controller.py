

import krpc
import time
import math
import numpy as np
import numpy.linalg as npl

import SC_solver as solver
import SC_params
from KSP_controller_utils import *

from threading import Thread

print('--------')

params = {}
with open('KSP_controller_params.txt', 'r', encoding='utf-8') as f:
    for line in f:
        pair = line.split('#')[0].split('=')
        if len(pair) == 2:
            key = pair[0].strip()
            value = eval(pair[1])
            params[key] = value

print('----------params------------')
for k in (params):
    print('  %s: \n%s' % (k, params[k]))
print('\n\ninitializing...')

deg2rad = math.pi / 180
rad2deg = 180 / math.pi
g0 = params['g0']

#连接krpc
print('connecting...')
conn = krpc.connect(name='SC_controller')
space_center = conn.space_center
vessel = space_center.active_vessel
flight = vessel.flight()
body = vessel.orbit.body

engine_gimbal = [m for m in vessel.parts.with_name('SSME')[0].modules if m.name == 'ModuleGimbal'][0]
# StarShip Main Engine(大雾)
engine_y = vessel.parts.with_name('SSME')[0].position(vessel.reference_frame)[1]

# starship flap
get_hinge = lambda tagname:[m for m in vessel.parts.with_tag(tagname)[0].modules if m.name=='ModuleRoboticServoHinge'][0]
h_fl = get_hinge('h_f_l')
h_fr = get_hinge('h_f_r')
h_rl = get_hinge('h_r_l')
h_rr = get_hinge('h_r_r')
#set 'Target Angle' 0~180
set_deploy = lambda h, deploy: h.set_field_float('Target Angle', math.asin(clamp(deploy, 0, 1)) * rad2deg + 5)
set_retract = lambda h:h.set_field_float('Target Angle', 0)

set_deploy(h_fl, 1)
set_deploy(h_fr, 1)
set_deploy(h_rl, 1)
set_deploy(h_rr, 1)

def combine_flaps(pitch_up, spin_right):
    pitch_up = clamp(pitch_up, -1, 1)
    #roll_right = clamp(roll_right, -1, 1)
    spin_right = clamp(spin_right, -1, 1)
    ctrl_fl = pitch_up
    ctrl_fr = pitch_up
    ctrl_rl = -pitch_up
    ctrl_rr = -pitch_up
    gap = max(0, 0.95 - max(ctrl_fl, ctrl_fr, ctrl_rl, ctrl_rr))
    ctrl_fl += gap
    ctrl_fr += gap
    ctrl_rl += gap
    ctrl_rr += gap
    ctrl_fl += -spin_right * 0.4
    ctrl_fr += spin_right * 0.4
    ctrl_rl += spin_right * 0.28
    ctrl_rr += -spin_right * 0.28
    set_deploy(h_fl, (ctrl_fl) / 2. + 0.5)
    set_deploy(h_fr, (ctrl_fr) / 2. + 0.5)
    set_deploy(h_rl, (ctrl_rl) / 2. + 0.5)
    set_deploy(h_rr, (ctrl_rr) / 2. + 0.5)

def retract_flaps():
    (h_fl).set_field_float('Target Angle', 5)
    (h_fr).set_field_float('Target Angle', 5)
    (h_rl).set_field_float('Target Angle', 5)
    (h_rr).set_field_float('Target Angle', 5)

#gimbal
#hinge_x = [m for m in vessel.parts.with_tag('hx')[0].modules if m.name=='ModuleRoboticServoHinge'][0]
#hinge_z = [m for m in vessel.parts.with_tag('hz')[0].modules if m.name=='ModuleRoboticServoHinge'][0]
#hinge_offset_y = hinge_x.part.position(vessel.reference_frame)[1]
#gimbalX = lambda angle:hinge_x.set_field_float('Target Angle', angle)
#gimbalY = lambda angle:hinge_z.set_field_float('Target Angle', angle)
#print(hinge_x.fields)
#www

delta_time = 0.01

#target
target_lat = params['target_lat'] * deg2rad
target_lon = params['target_lon'] * deg2rad
target_height = params['target_height']
target_axis = target_height + body.surface_height(target_lat * rad2deg, target_lon * rad2deg) + body.equatorial_radius
target_body_pos = np.array((math.cos(target_lon) * math.cos(target_lat), math.sin(target_lat), math.sin(target_lon) * math.cos(target_lat))) * target_axis


#limit
throttle_limit = params['throttle_limit']
throttle_limit_ctrl = params['throttle_limit_ctrl']
max_tilt = np.deg2rad(params['max_tilt'])
max_tilt_off = np.deg2rad(params['max_tilt_off'])

# krpc vessel对象生成vesselprofile
def get_vessel_profile(vessel):
    
    p = SC_params.VesselProfile()
    p.isp = vessel.specific_impulse
    p.g = vec(-g0, 0., 0.) # gravity
    p.m_dry = vessel.dry_mass
    p.gamma_gs = np.deg2rad(params['gamma_gs']) # glide slope
    p.theta_max = np.linspace(np.deg2rad(params['max_tilt']), np.deg2rad(10), SC_params.SuperParams().K) # tilt
    p.omega_max = np.deg2rad(params['max_omega']) # rotation vel
    p.delta_max = np.deg2rad(params['max_delta']) # gimbal
    p.T_min = vessel.available_thrust * throttle_limit[0]
    p.T_max = vessel.available_thrust * throttle_limit[1]
    p.r_T_B = vec(engine_y, 0., 0.) # thrust offset
    p.J_B_I = np.array(vessel.inertia_tensor).reshape((3, 3))
    p.airfric_k = params['airfric_k']
    p.time_guess = params['tf_guess']
    return p

# 根据krpc vessel对象预测落到给定高度时的vesselstate
def predict_vessel_state(vessel, est_height):
    # ref_target flight
    vel = vec(vessel.velocity(ref_target))
    pos = vec(vessel.position(ref_target))
    est_t = (pos[0] - est_height) / (-vel[0]) #匀速
    est_pos = pos + est_t * vel
    hdg_right = (flight.heading + 90) * deg2rad #右侧指向
    rot_axis = v3(0, math.cos(hdg_right), math.sin(hdg_right))
    rot_quat = quat(rot_axis, 90 * deg2rad)
    #qx, qy, qz, qw = vessel.rotation(ref_target) #xyzw转wxyz
    qx, qy, qz, qw = rot_quat #xyzw转wxyz
    
    state = SC_params.VesselState()
    state.mass = vessel.mass
    state.pos = est_pos
    state.vel = vel
    state.rotation = vec(qw, qx, qy, qz)
    #state.rotation = vec(1, 0, 0, 0)
    state.omega = vec(0, 0, 0)
    return state

def get_final_state(vessel, final_height):
    optimal_acc = vessel.available_thrust / vessel.mass * params['final_throttle'] - g0
    final_vel = math.sqrt(2 * optimal_acc * final_height)
    
    state = SC_params.VesselState()
    state.mass = vessel.mass
    state.pos = vec(final_height, 0, 0)
    state.vel = vec(-final_vel, 0, 0)
    state.rotation = vec(1, 0, 0, 0)
    state.omega = vec(0, 0, 0)
    return state

#fall attitude
ctrl_fall_pitch = PID()
ctrl_fall_pitch.kp = params['ctrl_fall_pitch.kp']
ctrl_fall_pitch.kd = params['ctrl_fall_pitch.kd']
#ctrl_fall_pitch.ki = params['ctrl_fall_pitch.ki'] # set that later
ctrl_fall_pitch.integral_limit = params['ctrl_fall_pitch.integral_limit']
ctrl_fall_yaw = PID()
ctrl_fall_yaw.kp = params['ctrl_fall_yaw.kp']
ctrl_fall_yaw.kd = params['ctrl_fall_yaw.kd']
ctrl_fall_distance = PID()
ctrl_fall_distance.kp = params['ctrl_fall_distance.kp']
ctrl_fall_distance.kd = params['ctrl_fall_distance.kd']
#rotation 
ctrl_x_rot = PID()
ctrl_x_rot.kp = params['ctrl_x_rot.kp']
ctrl_x_rot.kd = params['ctrl_x_rot.kd']
#ctrl_x_rot.redundancy = 0.1
ctrl_y_avel_kp = params['ctrl_y_avel_kp']
ctrl_z_rot = PID()
ctrl_z_rot.kp = params['ctrl_z_rot.kp']
ctrl_z_rot.kd = params['ctrl_z_rot.kd']
#ctrl_z_rot.redundancy = 0.1
#测量值
#torque = v3(3.66e+04, 5000, 3.66e+04)
#torque_k = v3(8.2e+04-3.66e+04, 0, 8.2e+04-3.66e+04)
#torque = v3(10300.000011920929, 10300.000011920929, 10300.000011920929)
#torque_k = v3(15183.20083618,  10772.2761631,   15183.24184418)
#print(vessel.available_torque)
# k
k_x = params['k_x']
k_v = params['k_v']

# final
final_throttle = params['final_throttle']
final_kp = params['final_kp']

# time init
game_delta_time = 0.02
game_prev_time = space_center.ut
start_time = time.time()

# references
print('creating target frame...')
ref_local = vessel.reference_frame 
ref_surface = vessel.surface_reference_frame #地面参考系
ref_body = body.reference_frame
ref_target_temp = space_center.ReferenceFrame.create_relative(ref_body, position=target_body_pos)
ref_target = space_center.ReferenceFrame.create_hybrid(ref_target_temp, rotation=ref_surface, velocity=ref_target_temp)

prev_vel = vec(vessel.velocity(ref_surface))
K = SC_params.SuperParams().K
solved_path = None
n_i = -1
error = vec(vessel.position(ref_target))
print('current error: %s' % error)
debug_lines = params['debug_lines']
if debug_lines:
    print('debug lines...')
    lines = [conn.drawing.add_line((0,0,0),(0,0,0), ref_target) for i in range(K-1)]
    directions = [conn.drawing.add_line((0,0,0), (1,0,0), ref_target) for i in range(K)]
    thrustvecs = [conn.drawing.add_line((0,0,0), (1,0,0), ref_target) for i in range(K)]
    target_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    target_line.color = (0,0,1)
    target2_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    target2_line.color = (0,0,1)
    head_line = conn.drawing.add_line((0,0,0),(1,0,0),ref_target)
    head_line.color = (0,1,1)
    for line in directions:
        line.color = (1,0,0)
    for line in thrustvecs:
        line.color = (1,0,1)

nav_mode = 'none'
frcount = 0

def update_lines(x, u):
    print('debug lines...')
    m_u = vessel.available_thrust
    for i in range(K-1):
        lines[i].start = x[1:4, i]
        lines[i].end = x[1:4, i+1]
    for i in range(K):
        mat = rotation_mat(x[7:11, i])
        directions[i].start = x[1:4, i]
        directions[i].end = x[1:4, i] + transform(vec(1, 0, 0), mat) * 5
        thrustvecs[i].start = x[1:4, i]
        thrustvecs[i].end = x[1:4, i] - transform(u[:, i], mat) / m_u * 10

def find_nearest_index(rk, vk, error):
    nearest_mag = npl.norm(rk[:, 0] - error)
    nearest_i = 0
    for i in range(x.shape[1]):
        mag = npl.norm(rk[:, i] - error) # + npl.norm(x[3:6, i] - v) * 0.2
        if mag < nearest_mag:
            nearest_mag = mag
            nearest_i = i
    v = vk[:, nearest_i]
    v_norm = npl.norm(v)
    v_dir = v / v_norm
    frac = clamp(np.dot(error - rk[:, nearest_i], v_dir) / (tf / K * v_norm), 0.5, -0.5)
    return nearest_i + frac

def sample_index(index, rk, vk, qk, uk):
    #if index >= N-1:
    if index >= K-1:
        return (rk[:, K-1], vk[:, K-1], qk[:, K-1], uk[:, K-1])
        #return (v3(0,0,0), v3(0,0,0), v3( 9.807,0,0))
    elif index <= 0:
        i = 0
        frac = index
    else:
        i = math.floor(index)
        frac = index - i
    r_i_s = lerp(rk[:, i], rk[:, i+1], frac)
    v_i_s = lerp(vk[:, i], vk[:, i+1], frac)
    q_i_s = lerp(qk[:, i], qk[:, i+1], frac)
    u_i_s = lerp(uk[:, i], uk[:, i+1], frac)
    if index < 0:
        u_i_s = uk[:, 1].copy()
        #print('u1234 ' + str(u[:, 0:4]))
        #print('u_i_s ' + str(u[:, 1]))
    return (r_i_s.copy(), v_i_s.copy(), q_i_s.copy(), u_i_s.copy())

def conic_clamp(target_a, min_mag, max_mag, max_tilt):
    a_mag = npl.norm(target_a)
    hor_dir = v3(0, target_a[1], target_a[2])
    hor_dir /= npl.norm(hor_dir)
    #target_direction = target_a / a_mag
    a_hor = npl.norm(target_a[1:3])
    a_ver = target_a[0]
    
    if (a_hor < min_mag * math.sin(max_tilt)):
        a_ver_min = math.sqrt(min_mag**2 - a_hor**2)
    else:
        a_ver_min = math.cos(max_tilt) * min_mag
    
    if (a_hor < max_mag * math.sin(max_tilt)):
        a_ver_max = math.sqrt(max_mag**2 - a_hor**2)
    else:
        a_ver_max = math.cos(max_tilt) * max_mag
    
    a_ver = clamp(a_ver, a_ver_max, a_ver_min)
    
    a_hor = min(a_hor, a_ver * math.tan(max_tilt))
    
    return hor_dir * a_hor + v3(a_ver, 0, 0)

def solve_path(vessel_profile, vessel_state, vessel_final_state):
    global solved_path, n_i
    print('----------vessel_profile(original)------------')
    for k in (vessel_profile.__dict__):
        print('-- %s: \n%s' % (k, vessel_profile.__dict__[k]))
    print('----------vessel_state(original)------------')
    for k in (vessel_state.__dict__):
        print('-- %s: \n%s' % (k, vessel_state.__dict__[k]))
    print('----------vessel_final_state(original)------------')
    for k in (vessel_final_state.__dict__):
        print('-- %s: \n%s' % (k, vessel_final_state.__dict__[k]))
        
    solver_options = SC_params.SolverOptions()
    solver_options.w_delta = lambda i:(1e-3 * (2 ** i))
    #solver_options.w_nu = 1e5
        
    print('---------solving----------')
    solved_path = solver.solve(vessel_profile, vessel_state, vessel_final_state, 
        solver_options=solver_options, use_c=True, verbose=True)
    if solved_path != None:
        (x, u, tf) = solved_path
        qw, qx, qy, qz = x[7:11, :]
        x[7:11, :] = vec(qx, qy, qz, qw) # wxyz转xyzw
            
        n_i = -100
        solved_path = (x, u, tf)
#        print('x slice')
#        print(x[:, 0:3])
#        print('u slice')
#        print(u[:, 0:3])
        print('---------solve done----------')
        if debug_lines:
            update_lines(x, u)
    else:
        print('---------solve error----------')

print('---------loop start-------------')
combine_flaps(0, 0)
while True:
    time.sleep(delta_time)
    real_time = time.time() - start_time
    ut = space_center.ut
    game_delta_time = ut - game_prev_time
    if game_delta_time < 0.01: #意味着游戏中还没有经过一个物理帧，所以不进行计算
        continue
        
# 取得一些之后要用的数据
    vessel_d = {}
    error = vessel_d['error'] = vec(vessel.position(ref_target)) # 目标系里的偏差
    avel = vessel_d['avel'] = vec(vessel.angular_velocity(ref_surface)) # 地面系下角速度（等于目标系角速度
    vel = vessel_d['vel'] = vec(vessel.velocity(ref_target)) # 地面速度
    
    rotation_local2srf = rotation_mat(vec(vessel.rotation(ref_surface))) # 机体系到地面系旋转矩阵
    rotation_srf2local = npl.inv(rotation_local2srf) # 地面系到机体系旋转矩阵
    moment_of_inertia_local = vec(vessel.moment_of_inertia) # 转动惯量
    mass = vessel_d['mass'] = vessel.mass
    max_thrust = vessel_d['max_thrust'] = vessel.available_thrust
    acceleration = vessel_d['acceleration'] = (vel - prev_vel) / game_delta_time
    
    #print(game_delta_time)
    
        
    if nav_mode == 'launch': #跳到一定高度
        balanced_thr = mass * g0 / max_thrust
        target_direction = v3(1, 0.02, 0) #偏北一点
        #print(target_direction)
        vessel.control.throttle = balanced_thr + (params['hop_vel'] - npl.norm(vel)) * 0.05
        #print(error[0])
        if (error[0] > params['hop_altitude']):
            nav_mode = 'transit'
            print('transit')
    elif nav_mode == 'transit': #减弱推力直到垂直速度为负
        balanced_thr = mass * g0 / max_thrust
        target_direction = v3(1, 0, 0)
        vessel.control.throttle = balanced_thr * 0.25
        if (vel[0] < -10):
            vessel.control.rcs = False
            vessel.control.pitch = -1
            time.sleep(1) #保持一秒
            vessel.control.pitch = 0
            vessel.control.throttle = 0
            nav_mode = 'fall'
            print('fall')
    elif nav_mode == 'fall': # 下落阶段 翼面控制姿态
        pitch_target = clamp(ctrl_fall_distance.update((math.sqrt(error[2]**2 + error[1]**2) - params['ctrl_fall_distance_target']) / 200., game_delta_time), -1, 1) * 15
        pitch_error = (flight.pitch - pitch_target) * deg2rad
        hdg_target = math.atan2(-error[2], -error[1]) * rad2deg
        hdg_error = norm_deg(flight.heading - hdg_target) * deg2rad
        #print(ctrl_fall_pitch.integral, math.sqrt(error[2]**2 + error[1]**2))
        if (abs(pitch_error) < 0.3):
            ctrl_fall_pitch.ki = params['ctrl_fall_pitch.ki']
        pitch_flap = ctrl_fall_pitch.update(pitch_error, game_delta_time)
        yaw_flap = ctrl_fall_yaw.update(hdg_error, game_delta_time)
        combine_flaps(pitch_flap, yaw_flap) #+0.1trim
        if error[0] < params['start_altitude']: #开始规划路径
            #frcount -= 1
            if frcount <= 0:
                frcount = 1
                vessel_profile = get_vessel_profile(vessel)
                vessel_state = predict_vessel_state(vessel, params['predict_altitude']) #计算要耗时，所以代入外推预测的未来状态
                vessel_final_state = get_final_state(vessel, params['final_height'])
                conn.krpc.paused = True  # ザ·ワールド
                #Thread(target=solve_path, args=(vessel_profile, vessel_state)).start()
                solve_path(vessel_profile, vessel_state, vessel_final_state)
                conn.krpc.paused = False
        if (error[0] <= params['predict_altitude'] and solved_path != None):
            n_i = 0
            vessel.control.sas = False
            vessel.control.rcs = True
            nav_mode = 'convex'
            retract_flaps()
            print('convex')
    elif nav_mode == 'convex': #沿路径
        (x, uk, tf) = solved_path
        mk = x[0, :] #mass
        rk = x[1:4, :] # position
        vk = x[4:7, :] # vel
        qk = x[7:11, :] # quaternion
        wk = x[11:14, :] # omega
            
        di = game_delta_time * K/tf
        n_i = clamp(find_nearest_index(rk, vk, error), n_i + di * 0.5, n_i + di * 1.5)
        #n_i = max(n_i - game_delta_time * 0.2 * K/tf, find_nearest_index(rk, vk, error)) #找到当前最近轨迹位置
        #print(game_delta_time)
        (r_i, v_i, q_i, u_i) = sample_index(n_i, rk, vk, qk, uk) #规划的位置速度
        (r_i_, v_i_, q_i_, u_i_) = sample_index(n_i + 0.4 * K/tf, rk, vk, qk, uk) #预测一小段时间以后
        q_i_mat = rotation_mat(q_i)
        q_i_mat_ = rotation_mat(q_i_)
        u_i = transform(u_i, q_i_mat)
        u_i_ = transform(u_i_, q_i_mat_)
        head_i = transform(vec(1, 0, 0), q_i_mat)
        head_i_ = transform(vec(1, 0, 0), q_i_mat_)
        #v_i_dir = v_i / npl.norm(v_i)
#        target_a = u_i_
#        target_v = v_i_
#        target_x = r_i # + np.dot((error - r_i), v_i_dir) * v_i_dir
#        #print(n_i, target_a, target_v, target_x)
#        target_a += (target_v - vel) * k_v + (target_x - error) * k_x
        
#        target_a = u_i + (v_i - vel) * k_v + (r_i - error) * k_x
#        target_a_ = u_i_ + (v_i_ - vel) * k_v + (r_i - error) * k_x
        target_a = npl.norm(u_i) / mass * head_i + (v_i - vel) * k_v + (r_i - error) * k_x
        target_a_ = npl.norm(u_i_) / mass * head_i_ + (v_i - vel) * k_v + (r_i - error) * k_x
        
        if debug_lines:
            target_line.start = error
            target_line.end = (r_i[0], r_i[1], r_i[2])
            target2_line.start = error
            target2_line.end = (r_i_[0], r_i_[1], r_i_[2])
        
        max_throttle_ctrl = throttle_limit_ctrl[1] * (max_thrust / mass)
        min_throttle_ctrl = throttle_limit_ctrl[0] * (max_thrust / mass)
        target_a = transform(target_a, q_i_mat.T)
        target_a = conic_clamp(target_a, min_throttle_ctrl, max_throttle_ctrl, max_tilt_off)
        target_a = transform(target_a, q_i_mat)
        target_a_ = transform(target_a_, q_i_mat_.T)
        target_a_ = conic_clamp(target_a_, min_throttle_ctrl, max_throttle_ctrl, max_tilt_off)
        target_a_ = transform(target_a_, q_i_mat_)
#        if n_i < 0 :
#            target_a = np.array([g0, 0, 0]) + u_i
        
        #target_direction = target_a_ / npl.norm(target_a_)
        #target_throttle = npl.norm(target_a) / (max_thrust / mass)
        target_direction = target_a_ / npl.norm(target_a_)
        target_throttle = npl.norm(target_a) / (max_thrust / mass)
        #print(target_a)
        if debug_lines:
            head_line.start = error
            head_line.end = error + target_direction * 8
        
        if n_i > 0:
            vessel.control.throttle = target_throttle
        
        if (K - n_i) * tf / K < 6:
            vessel.control.gear = True
        if npl.norm(error[1:3]) < params['final_radius'] and npl.norm(error[0]) < params['final_height']:
            vessel.control.gear = not vessel.control.gear
            engine_gimbal.set_field_float('Gimbal Limit', 30) # 防止震荡
            print('final')
            nav_mode = 'final'
    elif nav_mode == 'final':
        max_acc = throttle_limit_ctrl[1] * (max_thrust / mass) - g0
        max_acc_low = throttle_limit_ctrl[1] * final_throttle * (max_thrust / mass) - g0
        est_h = error[0] - vel[0]**2 / (2 * max_acc)
        est_h_low = error[0] - vel[0]**2 / (2 * max_acc_low)
        est_h_center = (est_h + est_h_low) / 2
        
        vessel.control.throttle = clamp(lerp(throttle_limit_ctrl[1] * final_throttle, throttle_limit_ctrl[1], -est_h_low / (est_h - est_h_low) * (1+final_kp)), throttle_limit_ctrl[1], throttle_limit_ctrl[0])
        
        error_hor = v3(0, error[1], error[2])
        vel_hor = v3(0, vel[1], vel[2])
        ctrl_hor = -error_hor * 0.03 - vel_hor * 0.06
        target_direction = ctrl_hor + v3(1, 0, 0)
        target_direction /= npl.norm(target_direction)
        target_direction = conic_clamp(target_direction, 1, 1, max_tilt)
    else:
        nav_mode = 'launch'
        if (max_thrust == 0):
            vessel.control.activate_next_stage()
        vessel.control.rcs = True
        vessel.control.gear = not vessel.control.gear
        continue
        #target_direction = -vel
        #vessel.control.throttle = 0
        
        
# 变换到机体坐标系计算姿态控制，以下xyz均指机体系
    if nav_mode in ['final', 'convex', 'launch', 'transit']:
        target_direction_local = transform(target_direction, rotation_srf2local) # 机体系的目标姿态的机体y轴指向
        avel_local = transform(avel, rotation_srf2local)  # 机体系角速度
        # 三个轴方向能提供的最大角加速度
        #authority_local = (torque + torque_k * vessel.control.throttle) / moment_of_inertia_local 
        #authority_local = np.abs(vec(vessel.available_torque[0])) / moment_of_inertia_local 
        #ctrl_x_rot.unified_authority = authority_local[0]
        #ctrl_z_rot.unified_authority = authority_local[2]
        # pid控制，roll直接消除角速度
        #control_pitch = -clamp(ctrl_x_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(1, 0, 0)), avel_local[0]), 1, -1)
        #control_yaw = -clamp(ctrl_z_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(0, 0, 1)), avel_local[2]), 1, -1)
        control_pitch = -clamp(ctrl_x_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(1, 0, 0)), game_delta_time), 1, -1)
        control_yaw = -clamp(ctrl_z_rot.update(angle_around_axis(target_direction_local, v3(0, 1, 0), v3(0, 0, 1)), game_delta_time), 1, -1)
        control_roll = clamp(avel_local[1] * ctrl_y_avel_kp, 1, -1)
        
        vessel.control.pitch = control_pitch
        vessel.control.yaw = control_yaw
        vessel.control.roll = control_roll
    
# 终止条件
    if nav_mode == 'final':
        if (npl.norm(error[1:3]) < 3 and npl.norm(error[0]) < 1 and npl.norm(vel[1:3]) < 0.3 and npl.norm(vel[0]) < 0.5 and npl.norm(avel) < 0.2) or (vel[0] > 0 and npl.norm(error[0]) < 1):
            print('exit')
            vessel.control.throttle = 0
            break
    
    prev_vel = vel
    game_prev_time = ut
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import pickle


def plot_X(X_in, U_in):

    ''' state variable
                   0   1   2   3   4   5   6   7   8   9  10  11  12  13
              x = [m, r0, r1, r2, v0, v1, v2, q0, q1, q2, q3, w0, w1, w2].T

        control variable, thrust
                  u = [T0, T1, T2] in body axes
    '''

    x = np.array(X_in)
    u = np.array(U_in)

    T = np.array([np.linalg.norm(uk) for uk in u])

    m = x[:,0]
    r = x[:,1:4]
    v = x[:,4:7]
    q = x[:,7:11]
    w = x[:,11:14]

    vm = np.array([np.linalg.norm(vk) for vk in v])

    qm = np.array([1. - 2.*(qk[2]**2 + qk[3]**2) for qk in q])
    qm = np.arccos(qm)
    qm = np.rad2deg(qm)

    dm = np.array([uk[0]/np.linalg.norm(uk) for uk in u])
    dm = np.arccos(dm)
    dm = np.rad2deg(dm)

    wm = np.array([wk for wk in w])
    wm = np.rad2deg(wm)

    plt.figure(1)

    plt.subplot(321)
    plt.plot(   r[:,1], r[:,0], 'k', r[:,1], r[:,0],'g.', label='Up/East')
    plt.quiver( r[:,1], r[:,0], u[:,1], u[:,0], label='Control Vectors',
                scale=50,color='b',headwidth =2,pivot='tip')
    plt.legend()

    plt.subplot(322)
    plt.plot(qm, 'k', qm, 'b.', label='Tilt Angle [deg]')
    plt.plot((90,)*len(qm))
    plt.legend()

    plt.subplot(323)
    plt.plot( T, 'k', T, 'c.', label='Thrust Magnitude')
    plt.legend()

    plt.subplot(324)
    plt.plot(wm, 'k', wm, 'b.', label='Angular Rate [deg/UT]')
    plt.plot((+60,)*len(wm))
    plt.legend()

    plt.subplot(325)
    plt.plot( m, 'k', m, 'm.', label='Mass')
    plt.legend()

    plt.subplot(326)
    plt.plot(dm,label='Thrust Pointing Angle [deg]')
    plt.plot((20,)*len(dm))
    plt.legend()
    plt.show()


def cIB(q):
    q0, q1, q2, q3 = q

    cIB_m = np.zeros((3,3))

    cIB_m[0,0] = 1-2*(q2**2 + q3**2)
    cIB_m[0,1] = 2*(q1*q2 + q0*q3)
    cIB_m[0,2] = 2*(q1*q3 - q0*q2)

    cIB_m[1,0] = 2*(q1*q2 - q0*q3)
    cIB_m[1,1] = 1-2*(q1**2 + q3**2)
    cIB_m[1,2] = 2*(q2*q3 + q0*q1)

    cIB_m[2,0] = 2*(q1*q3 + q0*q2)
    cIB_m[2,1] = 2*(q2*q3 - q0*q1)
    cIB_m[2,2] = 1-2*(q1**2 + q2**2)

    return cIB_m

def omega_mat(w):
    wx, wy, wz = w
    return np.mat([
        [0, -wx, -wy, -wz],
        [wx, 0, wz, -wy],
        [wy, -wz, 0, wx],
        [wz, wy, -wx, 0]
    ])

def rotation_mat(q):
    q = q / np.linalg.norm(q)
    #print(np.linalg.norm(q))
    (w, x, y, z) = q
    return np.mat([
    [1-2*y**2-2*z**2, 2*x*y+2*w*z, 2*x*z-2*w*y],
    [2*x*y-2*w*z, 1-2*x**2-2*z**2, 2*y*z+2*w*x],
    [2*x*z+2*w*y, 2*y*z-2*w*x, 1-2*x**2-2*y**2]
    ]).T

def transform(vec, mat):
    res = mat * np.mat(vec).T
    res = res.T
    #res = np.mat(vec) * mat
    return np.array((res[0,0], res[0,1], res[0,2]))

def plot_3D(X_in, U_in):

    x = np.array(X_in)
    u = np.array(U_in)

    r = x[:,1:4]
    v = x[:,4:7]
    q = x[:,7:11]
    omega = x[:, 11:14]
    #u[:, 1:3] *= -1
    rotations = [rotation_mat(q_) for q_ in q]
    u_rot = np.array([transform(u[i], rotations[i]) for i in range(len(rotations))])
    scale = 0.5
    localX = np.array([transform((scale, 0, 0), rotations[i]) for i in range(len(rotations))])
    localY = np.array([transform((0, scale, 0), rotations[i]) for i in range(len(rotations))])
    localZ = np.array([transform((0, 0, scale), rotations[i]) for i in range(len(rotations))])
    localomega = np.array([transform(omega[i], rotations[i]) for i in range(len(rotations))])
    
    acc = np.zeros((len(x), 3))
    acc[0] = 0
    dq_dt = np.zeros((len(x), 4))
    
    dt = 5.36411918 / 69.0
    #vel_integral = v[0]
    for i in range(1, len(u_rot) - 2):
        #vel_integral += dt * (u_rot[i] + u_rot[i+1]) / 2 / x[i, 0]
        #print(v[i+1], vel_integral)
        acc[i] = ((v[i+1] - v[i-1]) / 2. / dt + (1, 0, 0)) * 0.4
        dq_dt[i] = ((q[i+1] - q[i-1]) / 2. / dt)
        omegaw_q = 0.5 * (omega_mat(omega[i]) * np.mat(q[i]).T).T
#        if (i > 0 and i < 45):
#            print('dqdt\n', dq_dt[i], q[i], omega[i])
    
    u_rot *= 0.2
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    dord = (1, 2, 0)
    
    colors = []
    for i, _ in enumerate(r[:, 0]):

        vnorm = np.linalg.norm(v[i, :]) / 6.16

        colori = plt.cm.plasma(vnorm)
        colors.append(colori)

        ax.plot(
                        r[i : i+2, dord[0]],
                        r[i : i+2, dord[1]],
                        r[i : i+2, dord[2]],
                        color = colori,
                )


    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
                u_rot[:,dord[0]], u_rot[:,dord[1]], u_rot[:,dord[2]],

                pivot='tip',
                color=(1,0,1),
                arrow_length_ratio=0)
#    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
#                acc[:,dord[0]], acc[:,dord[1]], acc[:,dord[2]],
#                pivot='tip',
#                color=(1,0.4,0),
#                arrow_length_ratio=0)
    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
                localomega[:,dord[0]], localomega[:,dord[1]], localomega[:,dord[2]],
                color=(0,0,0),
                arrow_length_ratio=0)
                
    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
                localX[:,dord[0]], localX[:,dord[1]], localX[:,dord[2]],
                color=(1,0,0),
                arrow_length_ratio=0)
#    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
#                localY[:,dord[0]], localY[:,dord[1]], localY[:,dord[2]],
#                color=(0,1,0),
#                arrow_length_ratio=0)
#    ax.quiver(  r[:,dord[0]], r[:,dord[1]], r[:,dord[2]],
#                localZ[:,dord[0]], localZ[:,dord[1]], localZ[:,dord[2]],
#                color=(0,0,1),
#                arrow_length_ratio=0)
    
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([r[:,1].max()-r[:,1].min(), r[:,2].max()-r[:,2].min(), r[:,0].max()-r[:,0].min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(r[:,dord[0]].max()+r[:,dord[0]].min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(r[:,dord[1]].max()+r[:,dord[1]].min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(r[:,dord[2]].max()+r[:,dord[2]].min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')
    
    
    ax.set_aspect('equal')
    plt.show()

if __name__ == '__main__':
    X_in = pickle.load(open("X.p", "rb"))
    U_in = pickle.load(open("U.p", "rb"))
    plot_3D(X_in, U_in)
    plot_X(X_in, U_in)

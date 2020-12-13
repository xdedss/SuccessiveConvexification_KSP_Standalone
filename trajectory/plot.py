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

def plot_3D(X_in, U_in):

    x = np.array(X_in)
    u = np.array(U_in)

    r = x[:,1:4]
    v = x[:,4:7]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    colors = []
    for i, _ in enumerate(r[:, 0]):

        vnorm = np.linalg.norm(v[i, :]) / 6.16

        colori = plt.cm.plasma(vnorm)
        colors.append(colori)

        ax.plot(
                        r[i : i+2, 1],
                        r[i : i+2, 2],
                        r[i : i+2, 0],
                        color = colori,
                )


    ax.quiver(  r[:,1], r[:,2], r[:,0],
                u[:,1], u[:,2], u[:,0],

                pivot='tip',
                colors=colors,
                arrow_length_ratio=0)

    ax.set_aspect('equal')
    plt.show()

if __name__ == '__main__':
    X_in = pickle.load(open("X.p", "rb"))
    U_in = pickle.load(open("U.p", "rb"))
    plot_3D(X_in, U_in)
    plot_X(X_in, U_in)


import numpy as np
from time import time
from scipy.integrate import odeint
from scipy.special import gamma, airy
from numba import jit

y1_0 = 1.0 / 3**(2.0/3.0) / gamma(2.0/3.0)
y0_0 = -1.0 / 3**(1.0/3.0) / gamma(1.0/3.0)
y0 = [y0_0, y1_0]
t = np.arange(0, 4.0, 0.01)

def RHS(y, t):
    return [t*y[1],y[0]]

@jit(nopython=True)
def RHS2(y, t):
    y[0],y[1] = t*y[1],y[0]
    return y

start_time = time()
for i in range(1000):
    res = odeint(RHS, y0, t)
print('time=', time() - start_time)

start_time = time()
for i in range(1000):
    res = odeint(RHS2, y0, t)
print('time=', time() - start_time)

print(RHS)
print(RHS2)
print("-------------")
print(RHS2.inspect_types())
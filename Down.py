# stage 2: falling down
from sympy import *
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

m = 0.1
g = 9.8
k = 9.11*10**-3
v0 = 300.0
t0 = 0
t = 0
h0 = 0

def dydx2(v):
    return g - k/m * v**0.999

def rungeKutta2(t0, v0, t, h):
    # Count number of iterations using step size or step height h
    n = (int)((t - t0) / h)
    # Iterate for number of iterations
    v = v0
    x = 0
    for i in range(1, n + 1):
        "Apply Runge Kutta Formulas to find next value of v"
        k1 = h * dydx2(v)
        k2 = h * dydx2(v + 0.5 * k1)
        k3 = h * dydx2(v + 0.5 * k2)
        k4 = h * dydx2(v + k3)
        # Update next value of v
        v = v + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        # Update next value of x: dx/dt = v(t)
        x = x + v * h

    return v,x

def bisection2(a,b,N):
    a_n = a
    b_n = b
    for n in range(1, N + 1):
        m_n = (a_n + b_n) / 2
        f_m_n = np.array(rungeKutta2(t0, v0, m_n, h)[1])-1720.9
        if (rungeKutta2(t0, v0, a_n, h)[1]-1720.9) * f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif (rungeKutta2(t0, v0, b_n, h)[1]-1720.9) * f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif abs(f_m_n) <= 0.0001:
            print("Found approximate solution.")
            print(m_n)
            return m_n
    print((a_n + b_n) / 2)
    return (a_n + b_n) / 2


v1 = 0
h = 0.2
print(rungeKutta2(t0, v1, 20, h))
print(rungeKutta2(t0, v1, 30, h))
# h=0 must lie in the range of t = [20,30])

a = 20
b = 30
t_down = bisection2(a,b,1000)
print("the time when h=0 is", t_down)
t_up = 14.8
print("the total time used is", t_up + t_down)

def y(t,r):
    y, v = r
    dy_dt = v
    dv_dt = (m*g - k*v**0.999) / m
    return dy_dt, dv_dt

def solve(t):
    time = np.linspace(0,t)
    sol = solve_ivp(y_down, (time[0],time[-1]), (-Highest, 0), t_eval = time)
    return -sol.y[0][-1]

# plot the graph
t_fallingdown = np.linspace(0, t_down)
sol_fallingdown = solve_ivp(y, (t_fallingdown[0],t_fallingdown[-1]), (0,v0), t_eval = t_fallingdown)
height = sol_fallingdown.y[0]
velocity = sol_fallingdown.y[1]

plt.figure(figsize = (12,8))
plt.plot(t_fallingdown, height, label = 'Height')
plt.plot(t_fallingdown, velocity, label = 'Velocity')

plt.xlabel('time: s')
plt.title('Falling down stage: Height vs. Velocity')
plt.legend()
plt.show()


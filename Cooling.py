from sympy import *
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

A = 33.33 # environment
T0 = 91.11 # body
t = 120
h = 0.1
t0 = 0

def dydx(k, T):
    return k*(A-T**1.01)

# Finds value of y for a given x using step size h
# and initial value y0 at x0.
def rungeKutta(t0, T0, x, h, k):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x - t0) / h)
    # Iterate for number of iterations
    T = T0
    for i in range(1, n + 1):
        "Apply Runge Kutta Formulas to find next value of y"
        k1 = h * dydx(k,T)
        k2 = h * dydx(k,  T + 0.5 * k1)
        k3 = h * dydx(k, T + 0.5 * k2)
        k4 = h * dydx(k, T + k3)

        # Update next value of y
        T = T + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


    #print("x is ", x0,"y is ", y )
    return T

def bisection(a,b,N):
    a_n = a
    b_n = b
    for n in range(1, N + 1):
        m_n = (a_n + b_n) / 2
        f_m_n = rungeKutta(t0, T0, t, h, m_n)-66.66
        if (rungeKutta(t0, T0, t, h, a_n)-66.66) * f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif (rungeKutta(t0, T0, t, h, b_n)-66.66) * f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            print(m_n)
            return m_n
        else:
            print("Bisection method fails.")
            return None
    print((a_n + b_n) / 2)
    return (a_n + b_n) / 2

print(rungeKutta(t0, T0, t, h, 0))
print(rungeKutta(t0, T0, t, h, 0.2))

a = 0
b = 0.2

k = bisection(a, b, 1000)
print('k is', k)

print("by using RK4, temperature is", rungeKutta(t0, T0, t, h, k))

def ForwardEuler(T0, N, h):
    """Solve u’=f(u,t), u(0)=U0, with n steps until t=T."""
    t = np.zeros(N+1)
    T = np.zeros(N+1) # T[n] is the solution at time t[n]
    T[0] = T0
    t[0] = 0
    for n in range(N):
        if T[n] >= 98.88:
            print(t[n], T[n])
            return t[n]
        t[n + 1] = t[n] + h
        T[n + 1] = T[n] - k * (A - T[n] ** 1.01)
    return t[-1]

print('the time person needs to be cooled is', ForwardEuler(91.11, 1000, 1))
print('the person will die at 11:', 60-ForwardEuler(91.11, 1000, 1))

def ForwardEuler2(T0, N, h):
    """Solve u’=f(u,t), u(0)=U0, with n steps until t=T."""
    t = np.zeros(N+1)
    T = np.zeros(N+1) # T[n] is the solution at time t[n]
    T[0] = T0
    t[0] = 0
    for n in range(N):
        if T[n] <= 44.44:
            print(t[n], T[n])
            return t[n]
        t[n + 1] = t[n] + h
        T[n + 1] = T[n] + k * (A - T[n] ** 1.01)
    return t[-1]

print('the time needed to decrease to 44.44F is', ForwardEuler2(91.11, 10000, 1))
print('the temperature will reach 44.44F at 5:', (ForwardEuler2(91.11, 1000, 1)-5*60))

# plot the graph
temp = []
cooling_time = ForwardEuler(91.11, 1000, 1)
time = np.linspace(0, cooling_time)
for t in time:
    temp.append(rungeKutta(t0, T0, t, h, k))
plt.figure(figsize = (12,8))
plt.plot(time, temp, label = 'Temperature')

plt.title('Temperature vs. Time')
plt.xlabel('Time: s')
plt.ylabel('Temperature')
plt.legend()
plt.show()



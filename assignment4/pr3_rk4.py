import numpy as np
import matplotlib.pyplot as plt

# Initial conditions kept const for the rest of the problem
k = 1
m = 1
ti = 0
tf = 12 * np.pi
x0 = 1
v0 = 0

# function of x at RHS
def f(x):
    return -k * x / m

# Define the fourth-order Runge-Kutta method
# N is the number of time division
def runge_kutta(f,N):
    dt = (tf - ti) / N
    t = np.linspace(ti, tf, N)
    x = np.zeros(N)
    v = np.zeros(N)
    #initial conditions
    x[0] = x0 
    v[0] = v0
    for i in range(1, N):
        #the coefficients of rk 4 method
        # change in position is a linear function of v so the coeff are written accordingly
        # change in velocity depends on the positions
        k1x = dt * v[i-1]  
        k1v = dt * f(x[i-1])
        k2x = dt * (v[i-1] + 0.5 * k1v)
        k2v = dt * f(x[i-1] + 0.5 * k1x)
        k3x = dt * (v[i-1] + 0.5 * k2v)
        k3v = dt * f(x[i-1] + 0.5 * k2x)
        k4x = dt * (v[i-1] + k3v)
        k4v = dt * f(x[i-1] + k3x)

        x[i] = x[i-1] + (1/6) * (k1x + 2 * k2x + 2 * k3x + k4x) # change of position    
        #change of velocity depending on the position (k(n) are function of position)
        v[i] = v[i-1] + (1/6) * (k1v + 2 * k2v + 2 * k3v + k4v) 
    return t, x, v

# Define the analytical solution for initial energy 
def E_analytical(t):
    return 0.5 * (k*x0**2 + m*v0**2)

# Define a function to calculate the computed total energy
def E_computed(x, v):
    return 0.5 * (k*x**2 + m*v**2)

# Define the number of mesh points
N_values = [1000, 10000]

# Loop over different values of N
for N in N_values:
    t, x, v = runge_kutta(f,N)
    E_analytic = E_analytical(t) #const over the time
    E_comp = E_computed(x, v)
    delta_E = E_analytic - E_comp
    plt.plot(t, delta_E, label=f"N = {N}")

# Plot the results
plt.xlabel('Time (t)')
plt.ylabel(' $\Delta E$ (E_analytic - E_computed)')
plt.legend()
plt.title('$\Delta E$ vs Time plot')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
k = 1
m = 1
ti = 0
tf = 12 * np.pi
x0 = 1
v0 = 0

# Define the differential equation of motion
def f(t, x, v):
    return -k * x / m

# Define the fourth-order Runge-Kutta method
def runge_kutta(f, ti, tf, x0, v0, N):
    dt = (tf - ti) / N
    t = np.linspace(ti, tf, N)
    x = np.zeros(N)
    v = np.zeros(N)
    x[0] = x0
    v[0] = v0
    for i in range(1, N):
        k1x = dt * v[i-1]
        k1v = dt * f(t[i-1], x[i-1], v[i-1])
        k2x = dt * (v[i-1] + 0.5 * k1v)
        k2v = dt * f(t[i-1] + 0.5 * dt, x[i-1] + 0.5 * k1x, v[i-1] + 0.5 * k1v)
        k3x = dt * (v[i-1] + 0.5 * k2v)
        k3v = dt * f(t[i-1] + 0.5 * dt, x[i-1] + 0.5 * k2x, v[i-1] + 0.5 * k2v)
        k4x = dt * (v[i-1] + k3v)
        k4v = dt * f(t[i-1] + dt, x[i-1] + k3x, v[i-1] + k3v)
        x[i] = x[i-1] + (1/6) * (k1x + 2 * k2x + 2 * k3x + k4x)
        v[i] = v[i-1] + (1/6) * (k1v + 2 * k2v + 2 * k3v + k4v)
    return t, x, v

# Define the analytical solution for initial energy
def E_analytical(t):
    return 0.5 * (x0**2 + v0**2 * np.cos(np.sqrt(k/m) * t)**2)

# Define a function to calculate the computed total energy
def E_computed(x, v):
    return 0.5 * (x**2 + v**2)

# Define the number of mesh points
N_values = [1000, 10000]

# Loop over different values of N
for N in N_values:
    t, x, v = runge_kutta(f, ti, tf, x0, v0, N)
    E_analytic = E_analytical(t)
    E_comp = E_computed(x, v)
    delta_E = E_analytic - E_comp
    plt.plot(t, delta_E, label=f"N = {N}")

# Plot the results
plt.xlabel('Time (t)')
plt.ylabel('Delta E (E_analytic - E_computed)')
plt.legend()
plt.title('Energy Conservation in Runge-Kutta')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Constants
lamb = 0.3  # Decay constant (sec^-1)
dt = 1e-4       # Time step (sec)
Nruns = 1000    # Number of runs/events

# Initial number of nuclei
N0_list = [10, 100, 1000]

# Function to simulate radioactive decay
def radioactive_decay(N0, lamb, dt, Nruns):
    times = np.arange(0, Nruns*dt, dt)  # Time array
    Nt = np.zeros(len(times))           # Number of nuclei array
    Nt[0] = N0                          # Set initial number of nuclei
    for i in range(1, len(times)):
        dN = -lamb * Nt[i-1] * dt    # Decay rate
        Nt[i] = Nt[i-1] + dN            # Update number of nuclei
    return times, Nt

# Perform simulations for different initial N0 values
for N0 in N0_list:
    avg_Nt = np.zeros(len(radioactive_decay(N0, lambda_, dt, Nruns)[1]))
    for j in range(Nruns):
        Nt = radioactive_decay(N0, lambda_, dt, Nruns)
        avg_Nt += Nt
    avg_Nt /= Nruns

    # _Plot t vs. ln(Nt)
    plt.plot(radioactive_decay(N0, lambda_, dt, Nruns)[0], np.log(avg_Nt), label=f'N0 = {N0}')

# Plot exact result
t_exact = np.linspace(0, Nruns*dt, num=1000)
N_exact = N0 * np.exp(-lambda_ * t_exact)
plt.plot(t_exact, np.log(N_exact), 'r', label='Exact')
plt.xlabel('Time (s)')
plt.ylabel('ln(N(t))')
plt.legend()
plt.show()

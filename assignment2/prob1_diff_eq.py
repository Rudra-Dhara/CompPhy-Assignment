import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy import sparse

#Here, a, b, and c are the diagonal, main diagonal, and subdiagonal entries, d is lhs of the eqn
#this code written by taking the help of the internet
def tdma_solver(a, b, c, d):
    n = len(b)
     

    c_dash = np.zeros(n)
    d_dash = np.zeros(n)
    # Perform forward elemination
    c_dash[0] = c[0] / b[0]
    for i in range(1, n):
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i-1])
        d_dash[i] = (d[i] - a[i] * d_dash[i-1]) / (b[i] - a[i] * c_dash[i-1])
    
    # Perform backward elimination
    x = np.zeros(n)
    x[n-1] = d_dash[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_dash[i] - c_dash[i] * x[i+1]
    
    return x





# Define the Poisson equation function
def poi_fn(x):
    return -(3*x + x**2)*np.exp(x)

# Define the exact solution
def exact(x):
    return x*(1-x)*np.exp(x)

# Define the grid and step size h
n_values = [10, 100, 1000, 10**4,10**6]
max_errors = []

for n in n_values:
    h = 1/(n+1)
    x = np.linspace(h, 1-h, n) #x values with different interval 

    #the 
    main_diag = -2*np.ones(n)
    off_diag = np.ones(n)   

    # b vector
    b = h**2*poi_fn(x)

    # Apply the boundary conditions
    b[0] = 0
    b[-1] = 0

    #solving u
    u = tdma_solver(off_diag, main_diag, off_diag, b)

    # Calculate the maximum relative error
    max_error = np.max(np.abs((u-exact(x))/exact(x)))
    max_errors.append(max_error)

    # Plot the numerical solution and the exact solution
    plt.plot(x, u, label='Numerical solution for {}'.format(n))
    
    
# Tabulate the results
plt.plot(x, exact(x), label='Exact solution')
plt.legend()
plt.title('comparision b/w exact solution')
plt.show()

print('h\tmax error')
for n, max_error in zip(n_values, max_errors):
    h = 1/(n+1)
    print(f'{h:.1e}\t{max_error:.8f}') 
    
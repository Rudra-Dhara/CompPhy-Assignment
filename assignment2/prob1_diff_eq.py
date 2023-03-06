import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy import sparse

# Define the Poisson equation function
def poi_fn(x):
    return -(3*x + x**2)*np.exp(x)

# Define the exact solution
def exact(x):
    return x*(1-x)*np.exp(x)

# Define the grid and step size h
n_values = [10, 100, 1000, 10000,100000, 1000000]
max_errors = []

for n in n_values:
    h = 1/(n+1)
    x = np.linspace(h, 1-h, n) #x values with different interval 

    
    main_diag = -2*np.ones(n)
    off_diag = np.ones(n)
    data = [off_diag, main_diag, off_diag]
    # Tridiagonal matrix
    A = sparse.spdiags(data,[-1,0,1],n,n,format='csc')
    

    # b vector of right side
    b = h**2*poi_fn(x)

    # Apply the boundary conditions
    b[0] = 0
    b[-1] = 0

    # Solve the linear system using scipy.linalg.solve
    u = spsolve(A, b)

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

import numpy as np

# Define the integrand function
def integrand(x):
    return np.exp(-x**2)

# Define the probability distribution function for importance sampling
def p(x, A):
    return A * np.exp(-x)

# Define the limits of integration
a = 0
b = 1

# Define the values of N for Monte Carlo sampling
N_values = [10**3, 10**4, 10**5]

# Loop over different values of N
for N in N_values:
    # (a) Brute force Monte Carlo
    x_brute = np.random.uniform(a, b, N) # Generate random numbers from uniform distribution
    integral_brute = np.mean(integrand(x_brute)) * (b - a) # Estimate the integral using the mean of function evaluations
    
    # (b) Importance sampling with p(x) = Ae^-x
    A = 1 / (1 - np.exp(-b)) # Calculate the value of A to normalize p(x)
    x_importance = np.random.exponential(scale=1/A, size=N) # Generate random numbers from exponential distribution with scale 1/A
    integral_importance = np.mean(integrand(x_importance) / p(x_importance, A)) * (b - a) # Estimate the integral using importance sampling
    
    # Print the results
    print(f"For N = {N}:")
    print(f"Brute force Monte Carlo: {integral_brute:.6f}")
    print(f"Importance sampling: {integral_importance:.6f}")
    print("-----")

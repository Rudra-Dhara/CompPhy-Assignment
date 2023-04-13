import numpy as np

# integrand function
def integrand(x):
    return np.exp(-x**2)

# PDF for importance sampling
def p(x, A):
    return A * np.exp(-x)

# integration limit
a = 0
b = 1

# values of N for Monte Carlo sampling
N_values = [10**3, 10**4, 10**5]

# loop over different values of N
for N in N_values:
    # (a) brute force Monte Carlo
    x_brute = np.random.rand(N) # Generate random numbers from uniform distribution
    integral_brute = np.mean(integrand(a+ x_brute*(b-a))*(b-a))  # Estimate the integral using the mean of function evaluations
    brute_var = np.var(integrand(a+ x_brute*(b-a))*(b-a))

    # (b) importance sampling with p(x) = Ae**-x
    # the transformed variable y(x) is = - ln (1 - x/A)
    A = 1 / (np.exp(-a) - np.exp(-b)) # Calculate the value of A to normalize p(x)
    x_imp = -np.log(np.ones_like(x_brute)*np.exp(-a) - x_brute/A) # change of variavle
    integral_imp = np.mean(integrand(x_imp) / p(x_imp, A)) # Estimate the integral using importance sampling
    imp_var= np.var(integrand(x_imp) / p(x_imp, A))
    # Print the results
    print(f"For N = {N}:")
    print(f"Brute force Monte Carlo: {integral_brute:.6f} with sd ={brute_var**0.5:6f}")
    print(f"Importance sampling: {integral_imp:.6f} with sd={imp_var**0.5:6f}")
    print("-----")

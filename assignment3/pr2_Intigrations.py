import numpy as np

#integrand function
def f(x):
    return 1/(2+x**2)

# trapezoidal rule
def trapezoidal(f, a, b, n):
    x = np.linspace(a, b, n+1)
    dx = (b-a)/n
    I = dx/2 * (f(a) + 2*np.sum(f(x[1:n])) + f(b))
    return I

#  Simpson's rule
def simpson(f, a, b, n):
    x = np.linspace(a, b, n+1)
    dx = (b-a)/n
    I = dx/3 * (f(a) + 4*np.sum(f(x[1:n:2])) + 2*np.sum(f(x[2:n-1:2])) + f(b))
    return I

#  Gauss-Legendre quadrature
def gauss_legendre(f, a, b, n):
    x, w = np.polynomial.legendre.leggauss(n+1)  #x is the point and w is the weight list of the function
    xp = (b-a)/2*x + (b+a)/2
    wp = (b-a)/2*w
    I = np.sum(wp*f(xp))
    return I



# Define the limits of integration
a, b = 0, 3

# mesh points
N = [10, 40, 100, 1000]

# integrals using the trapezoidal rule
for n in N:
    I = trapezoidal(f, a, b, n)
    print(f"Trapezoidal rule with {n} mesh points: I = {I:.6f}")
print('\n')

#using Simpson's rule
for n in N:
    I = simpson(f, a, b, n)
    print(f"Simpson's rule with {n} mesh points: I = {I:.6f}")
print('\n')

# Gauss-Legendre quadrature
for n in N:
    I = gauss_legendre(f, a, b, n)
    print(f"Gauss-Legendre quadrature with {n} mesh points: I = {I:.6f}")

import int_fn as int


#integrand function
def f(x):
    return 1/(2+x**2)
# Define the limits of integration
a, b = 0, 3

# mesh points
N = [10, 40, 100, 1000]

# integrals using the trapezoidal rule
for n in N:
    I = int.trapezoidal(f, a, b, n)
    print(f"Trapezoidal rule with {n} mesh points: I = {I:.6f}")
print('\n')

#using Simpson's rule
for n in N:
    I = int.simpson(f, a, b, n)
    print(f"Simpson's rule with {n} mesh points: I = {I:.6f}")
print('\n')

# Gauss-Legendre quadrature
for n in N:
    I = int.gauss_legendre(f, a, b, n)
    print(f"Gauss-Legendre quadrature with {n} mesh points: I = {I:.6f}")


print('\n\n')
# Define the integrand function
def integrand(x, y):
    return 1 / (1 - x*y)

# Define the limits of integration
x_min, x_max = 0, 0.9999        #as the integral is blowing up at x=1 we put the limit si
y_min, y_max = 0, 0.9999
#double integral's calculation
for n in N:
    I = int.trapz2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"Trapezoidal 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')

for n in N:
    I = int.simpson_2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"SImpson 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')

for n in N:
    I = int.gauss_legendre_2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"Gauss 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')    
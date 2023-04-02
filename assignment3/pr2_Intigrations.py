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






# Define the integrand function
def integrand(x, y):
    return 1 / (1 - x*y)

# Define the limits of integration
x_min, x_max = 0, 0.9999        #as the integral is blowing up at x=1 we put the limit si
y_min, y_max = 0, 0.9999

def trapz2d(f, a, b, c, d, nx=1000, ny=1000):
    """Approximate the integral of f(x,y) over the domain [a,b] x [c,d] using
    the composite trapezoidal rule with nx panels in x and ny panels in y."""
    
    # create the x and y mesh points
    x = np.linspace(a, b, nx+1)
    y = np.linspace(c, d, ny+1)
    
    # compute the step sizes in the x and y directions
    hx = (b - a) / nx
    hy = (d - c) / ny
    
    # initialize the sum
    integral = 0.0
    
    # loop over the panels in the x direction
    for i in range(nx):
        # loop over the panels in the y direction
        for j in range(ny):
            # evaluate the function at the four corners of the panel
            f00 = f(x[i], y[j])
            f01 = f(x[i], y[j+1])
            f10 = f(x[i+1], y[j])
            f11 = f(x[i+1], y[j+1])
            
            # compute the integral over this panel using the trapezoidal rule
            panel_integral = (hx*hy/4) * (f00 + f01 + f10 + f11)
            
            # add the panel integral to the total integral
            integral += panel_integral
    
    return integral

print(trapz2d(integrand,x_min,x_max,y_min,y_max))
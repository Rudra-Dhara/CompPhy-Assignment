import numpy as np
import math

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
def simpson(f,x0,xn,n,*args):
    # calculating step size
    h = (xn - x0) / n
    
    # Finding sum 
    integration = f(x0,*args) + f(xn,*args)
    
    for i in range(1,n):
        k = x0 + i*h
        
        if i%2 == 0:
            integration = integration + 2 * f(k,*args)
        else:
            integration = integration + 4 * f(k,*args)
    
    # Finding final integration value
    integration = integration * h/3
    
    return integration

#  Gauss-Legendre quadrature
def gauss_legendre(f, a, b, n,*args):
    x = [0.0] * n     # Roots of Legendre polynomial
    w = [0.0] * n     # Weights for each root
    xm = 0.5 * (b + a)   # Midpoint of interval [a,b]
    xr = 0.5 * (b - a)   # Half-length of interval [a,b]
    sum = 0.0     # Running sum of the integrand
    for i in range(n):
        m = (n + 1) // 2   # Index of midpoint
        # Compute Legendre polynomial using recurrence relation
        z = math.cos(math.pi * (i + 0.75) / (n + 0.5))   # Root of Legendre polynomial
        p1 = 1.0   # Legendre polynomial for i=0
        p2 = 0.0   # Legendre polynomial for i=-1 (not used)
        for j in range(1, n + 1):
            p3 = p2   # Legendre polynomial for i=j-2 (not used)
            p2 = p1   # Legendre polynomial for i=j-1
            # Recurrence relation for Legendre polynomial
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
        # Compute derivative of Legendre polynomial
        pp = n * (z * p1 - p2) / (z * z - 1.0)
        # Compute root and weight
        x[i] = xm - xr * z   # Root in interval [a,b]
        w[i] = 2.0 * xr / ((1.0 - z * z) * pp * pp)   # Weight for root
        # Update running sum
        sum += w[i] * f(x[i],*args)   # Weighted value of integrand at root
    return sum



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


print('\n\n')



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


def simpson_2d(f, a, b, c, d, Nx, Ny):
    hx = (b-a)/Nx
    hy = (d-c)/Ny
    
    x = np.linspace(a, b, Nx+1)
    y = np.linspace(c, d, Ny+1)
    
    
    integral = 0
    I=simpson(f,a,b,Nx,y[0])+simpson(f,a,b,Nx,y[-1])

    for i in range(1,Ny):
            if i%2==0:
                I += 2*simpson(f,a,b,Nx,y[i])
            else:
                I+= 4*simpson(f,a,b,Nx,y[i])
    I=I*hy/3

    return I

def gauss_legendre_2d(f, a, b, c, d, Nx, Ny):
    y, w = np.polynomial.legendre.leggauss(Ny+1)  #x is the point and w is the weight list of the function
    yp = (d-c)/2*y + (d+c)/2
    wp = (d-c)/2*w
    I = 0
    for i in range(Ny+1):
        I += wp[i]*gauss_legendre(f,a,b,Nx,yp[i])
    return I

for n in N:
    I = trapz2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"Trapezoidal 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')

for n in N:
    I = simpson_2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"SImpson 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')

for n in N:
    I = gauss_legendre_2d(integrand,x_min,x_max,y_min,y_max,n,n)
    print(f"Gauss 2d rule with {n} mesh points: I = {I:.6f}")
print('\n')    
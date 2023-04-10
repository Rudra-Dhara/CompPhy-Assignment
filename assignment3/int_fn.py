import numpy as np
import math

#this function is a general function to do intigral in theree different method
#you can do double integral using 2d methods 


# in the up coming function 
# -------------------------
# f - is the integrand function
# a - lower limit of x integral
# b - upper limit of x integral
# n - number of division points

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
# here the approximation for for root of legendre function is used and weight from derivative formula
def gauss_legendre(f, a, b, n,*args):
    x = [0] * n     # Roots of Legendre polynomial
    w = [0] * n     # Weights for each root
    xm = 0.5 * (b + a)   # Midpoint of interval [a,b]
    xr = 0.5 * (b - a)   # Half-length of interval [a,b]
    sum = 0     
    for i in range(n):
        # Compute Legendre polynomial using recurrence relation
        z = math.cos(math.pi * (i + 0.75) / (n + 0.5))   # Root of Legendre polynomial
        p1 = 1   # Legendre polynomial l=0
        p2 = 0   # Legendre polynomial l=-1

        #calculating the  n th odrer legendre
        for j in range(1, n + 1):
            p3 = p2   # Legendre polynomial for i=j-2 
            p2 = p1   # Legendre polynomial for i=j-1
            # Recurrence relation for Legendre polynomial 
            p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
        # Compute derivative of Legendre polynomial
        pp = n * (z * p1 - p2) / (z**2 - 1)       #taken from the book
        # Compute root and weight
        x[i] = xm - xr * z   # Root in interval [a,b]
        w[i] = 2 * xr / ((1 - z * z) * pp**2)   # Weight for root

        sum += w[i] * f(x[i],*args)   # Weighted value of integrand at root
    return sum


#defining function for double integrals
# c - lower limit of y integral
# d - upper limit of y integral
# nx and ny - are number of division along x and y 
# other arguments are same as above
def trapz2d(f, a, b, c, d, nx, ny):
    
    # create the x and y mesh points
    x = np.linspace(a, b, nx+1)
    y = np.linspace(c, d, ny+1)
    
    # step sizes in the x and y 
    hx = (b - a) / nx
    hy = (d - c) / ny
    
    integral = 0.0
    
    for i in range(nx):
        
        for j in range(ny):
            # function at the four corners of the trapizoid
            f00 = f(x[i], y[j])
            f01 = f(x[i], y[j+1])
            f10 = f(x[i+1], y[j])
            f11 = f(x[i+1], y[j+1])
            
            # volume of the box
            box_integral = (hx*hy/4) * (f00 + f01 + f10 + f11)
            
            integral += box_integral
    
    return integral

#here we considered the whole integral along strips along x direction 
# then we evaluate the total integral of strips along y direction
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
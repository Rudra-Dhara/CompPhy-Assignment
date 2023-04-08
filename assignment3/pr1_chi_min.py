import numpy as np
from scipy.optimize import minimize
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 

def gauss_legendre(f, a, b, n,*args):
    x = [0.0] * n     # Roots of Legendre polynomial
    w = [0.0] * n     # Weights for each root
    xm = 0.5 * (b + a)   # Midpoint of interval [a,b]
    xr = 0.5 * (b - a)   # Half-length of interval [a,b]
    sum = 0     
    for i in range(n):
        # Compute Legendre polynomial using recurrence relation
        z = math.cos(math.pi * (i + 0.75) / (n + 0.5))   # Root of Legendre polynomial
        p1 = 1.0   # Legendre polynomial l=0
        p2 = 0.0   # Legendre polynomial l=-1

        #calculating the  n th odrer legendre
        for j in range(1, n + 1):
            p3 = p2   # Legendre polynomial for i=j-2 
            p2 = p1   # Legendre polynomial for i=j-1
            # Recurrence relation for Legendre polynomial 
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
        # Compute derivative of Legendre polynomial
        pp = n * (z * p1 - p2) / (z**2 - 1.0)       #taken from the book
        # Compute root and weight
        x[i] = xm - xr * z   # Root in interval [a,b]
        w[i] = 2.0 * xr / ((1.0 - z * z) * pp**2)   # Weight for root

        sum += w[i] * f(x[i],*args)   # Weighted value of integrand at root
    return sum

# Define the density of states function
def rho(m, A, TH, m0=0.5):
    return A * (m**2 + m0**2)**(5/4) * np.exp(m/TH)

# Define the theoretical number of states function
def Ntheory(A, m0, TH, mmax):
    return gauss_legendre(lambda m: rho(m, A, m0, TH), 0, mmax, n=10)

# Define the step function
def step_fn(m, m_i):
    if m >= m_i:
        return 1
    else:
        return 0
# declaring the value of g as global variable

# Define the experimental number of states function for mesons
def Nexp_mesons(m,mes_list):
    g=3 # for meson
    sum = 0
    for i in mes_list:
        sum+= g*step_fn(m,i)

    return sum

# Define the experimental number of states function for baryons
def Nexp_baryons(m, bar_list):
    g= 4 #for baryon
    sum = 0
    for i in bar_list:
        sum+= g*step_fn(m,i)

    return sum

#defining chi**2 function as all sigma = 1 we eliminate it from the function
def chi2(Nexp,exp_list,A,TH,m0=0.5):
    sum=0
    for m_i in exp_list:
        sum+= (np.log10(Nexp(m_i,exp_list)) - np.log10(Ntheory(A, m0, TH, m_i)))**2
    
    return sum


# Define the experimental masses and degeneracies
mes_mass = np.array([0.782, 1.170, 1.282, 1.420, 1.512])
bar_mass = np.array([0.938, 1.440, 1.535, 1.650, 1.710])
sigma=np.ones(5)



# Define the chi-square function with fixed parameters
def chi2_fixed(params):
    A, TH = params
    bar_mass = np.array([0.938, 1.440, 1.535, 1.650, 1.710])
    return chi2(Nexp_baryons, bar_mass, A, TH, 0.5)







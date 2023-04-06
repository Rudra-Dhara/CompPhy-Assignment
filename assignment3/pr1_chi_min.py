import numpy as np
from pr2_Intigration.py import gauss_legendre

# Define the density of states function
def rho(m, A, m0=0.5, TH):
    return A * (m**2 + m0**2)**(5/4) * np.exp(m/TH)

# Define the theoretical number of states function
def Ntheory(A, m0, TH, mmax):
    return gauss_legendre(lambda m: rho(m, A, m0, TH), m0, mmax, n=10)

# Define the step function
def step(x):
    return 1.0 if x >= 0 else 0.0
# declaring the value of g as global variable

# Define the experimental number of states function for mesons
def Nexp_mesons(m,mes_list):
    g=3 # for meson
    sum = 0
    for i in mes_list:
        sum+= g*(m - i)

    return sum

# Define the experimental number of states function for baryons
def Nexp_baryons(m, bar_list):
    g= 4 #for baryon
    sum = 0
    for i in bar_list:
        sum+= g*(m - i)

    return sum

def chi2(Nexp,exp_list,m,sigma_list,A,TH,mmax,m0=0.5):
    sum=0
    for i in sigma_list:
        sum+= (np.log10(Nexp(m,exp_list)) - np.log10(Ntheory(A, m0, TH, mmax)))/i
    
    return sum


# Define the experimental masses and degeneracies
mes_mass = np.array([0.782, 1.170, 1.282, 1.420, 1.512])
bar_mass = np.array([0.938, 1.440, 1.535, 1.650, 1.710])
sigma=np.ones(5)


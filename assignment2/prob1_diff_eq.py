import numpy as np
import matplotlib.pyplot as plt

#Here, a, b, and c are the diagonal, main diagonal, and subdiagonal entries, d is lhs of the eqn
def tdma_solver(a, b, c, d):
    n = len(b)
     

    c_dash = np.zeros(n)
    d_dash = np.zeros(n)
    # Perform forward elemination
    c_dash[0] = c[0] / b[0]
    for i in range(1, n):
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i-1])
        d_dash[i] = (d[i] - a[i] * d_dash[i-1]) / (b[i] - a[i] * c_dash[i-1])
    
    # Perform backward elimination
    x = np.zeros(n)
    x[n-1] = d_dash[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_dash[i] - c_dash[i] * x[i+1]
    
    return x


# Define the Poisson equation function
def poi_fn(x):
    return -(3*x + x**2)*np.exp(x)

# Define the exact solution
def exact(x):
    return x*(1-x)*np.exp(x)

# Define the grid and step size h
n_values = [10, 100, 1000, 10**4,10**5,10**6]


#part a
print('With this program we solve')
list_of_errlst=[] # the lists to store the different values of u_list for different value
list_of_xlst=[]
for n in n_values:
    h = 1/(n)
    x = np.linspace(0,1, n+1) #x values with different interval 
    #the tridiagnolan matrix diagonals
    main_diag = -2*np.ones(n-1)
    off_diag = np.ones(n-1)   

    # b vector
    b = h**2*poi_fn(x)
    

    #solving u
    u = tdma_solver(off_diag, main_diag, off_diag, b)
    u= list(u)
    
    # relative error list and x (storing the different lists due to different n within the loop)
    u_err= []
    for i in range(len(u)):
        u_err.append(np.log10(abs((u[i]-exact(x[i+1])/exact(x[i+1])))))# plus one added as the boundary value of soln is not included here
                     
    list_of_errlst.append(u_err)
    list_of_xlst.append(x)

    #inserting the boundary values of u
    u.insert(0,0)
    u.append(0)

    # Plot the numerical solution and the exact solution
    plt.plot(x, u, label='Numerical solution for {}'.format(n))
    
    

# Part c
plt.plot(x, exact(x), label='Exact solution')
plt.legend()
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('comparision b/w exact solution')
plt.grid()
plt.show()


#Part d 

#finding the maximum relative error
max_err=[]
h_list=[]
for i in range(len(n_values)):
    h_list.append(np.log10(1/n_values[i]))
    max_err.append(max(list_of_errlst[i]))

plt.title('The max relative error(in log scale) vs x for different h in log scale')
plt.plot(h_list,max_err,marker='.',ls=':')
plt.xlabel('$log_{10}(h)$')
plt.ylabel('$log_{10}(rel err)$')
plt.legend()
plt.grid()
plt.show()


    



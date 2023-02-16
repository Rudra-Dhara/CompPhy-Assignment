import matplotlib.pyplot as plt
import numpy as np
import math as m

#finding the bessel function using up method
def j0(x):
    return np.sin(x)/x
def j1(x):
    return (np.sin(x)- x*np.cos(x))/x**2


def j_up(x,l):
    if l>1:
        return (2*l-1)*j_up(x,l-1)/x - j_up(x,l-2)
    elif l==1:    
        return j1(x)
    else:
        return j0(x)

#defining the function for down method
#initial function with garbage value (using asymptotic approximation)
def ini_j100(x,arr=np.array([])):
    return arr
    
def ini_j99(x,arr=np.array([])):
    return arr
    

#using for loop instead of recursion so the program will be fast
def ini_j_down(x,l,arr1=np.array([]),arr2=np.array([])):
    arr1=ini_j99(x,arr1)
    arr2= ini_j100(x,arr2)
    for i in range(100,l,-1):
        temp=arr1
        arr1= (2*i+1)*arr1/x - arr2
        arr2=temp
    return arr1
    
#main program

#different value of x
x=np.arange(0.1,50.2,0.2,np.longdouble)
i=10 # value of l

#for down method initial guess for large l spherical bessel function
j99_list=np.ones(len(x))
j100_list=np.ones(len(x))

#normalizing the guess of large bessel function
j99_list= (j1(x)/ini_j_down(x,1,j99_list,j100_list))
j100_list=j99_list
  


print("x vs j{}(x)up and j{}(x)down is given below(respectively):".format(i,i))
for k in range(len(x)):
    print("{:2f}   {:6e}   {:6e}".format(x[k],j_up(x[k],i),ini_j_down(x,i,j99_list,j100_list)[k]))


#error calculation
err_arr=abs(j_up(x,i)-ini_j_down(x,i,j99_list,j100_list))/(abs(j_up(x,i))+abs(ini_j_down(x,i,j99_list,j100_list)))

print("\n\n\n The value of x vs relative error is given below respectively")
for k in range(len(x)):
    print("{:.2f}   {:6e}".format(x[k],err_arr[k]))

#log-log plotting the bessel function in different method
plt.title('log-log plot of bessel functions ')
plt.plot(np.log(ini_j_down(x,i,j99_list,j100_list)),np.log(x),label="$j_{}(x)$ in DOWN method".format({i}))
plt.plot(np.log(j_up(x,i)),np.log(x),label="$j_{}(x)$ in UP method".format({i}))
plt.ylabel('log(x)')
plt.xlabel('value of $log(j_{}(x))$'.format({i}))
plt.legend()
plt.grid()
plt.show()

plt.title("log log plot of relative errors")
plt.plot(np.log(x),np.log(err_arr),label='relative error')
plt.xlabel('log(x)')
plt.ylabel('log($\epsilon$)')
plt.legend()
plt.grid()
plt.show()
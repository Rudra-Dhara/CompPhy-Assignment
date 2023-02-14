import matplotlib.pyplot as plt
import numpy as np
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
    
x=np.arange(0.2,50.2,0.2,np.longdouble)
i=10


print("x vs j{}(x) is given below:".format(i))
for k in x:
    print("{}   {}".format(round(k,2),j_up(k,i)))


plt.plot(x,j_up(x,i),label="j{}(x)".format(i))
plt.xlabel('x')
plt.ylabel('value of j{}(x)'.format(i))
plt.legend()
plt.grid()
plt.show()
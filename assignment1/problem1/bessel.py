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
    
x=np.linspace(0,50,250)
i=0
plt.plot(x,j_up(x,i),label="{}".format(i))
plt.legend()
plt.show()




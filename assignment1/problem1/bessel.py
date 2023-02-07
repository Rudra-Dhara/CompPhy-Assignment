import matplotlib.pyplot as plt
import numpy as np
def j0(x):
    return np.sin(x)/x
def j1(x):
    return (np.sin(x)- x*np.cos(x))/x**2
def dF(f,x,*arg):
    h=1
    return (f(x-2*h,*arg)-8*f(x-h,*arg)+8*f(x+h,*arg)-f(x+2*h,*arg))/(12*h)

def j_up(x,l):
    if l>1:
        return (l-1)*j_up(x,l-1)/x - dF(j_up,x,l-1)
    elif l==1:    
        return j1(x)
    else:
        return j0(x)
    
x=np.linspace(0.1,50,100)
for i in x:
    print(i,"   ",j_up(i,10))

import math as m
#defining the function for bisection method
def bisec_rt(f,a,b,x_tol):
    h=x_tol
    if f(a)*f(b)>0:
        print("please enter a valid input")
    elif f(a)*f(b)<0:
        c= (a+b)/2
        count=0
        while abs(a-b)>=h:
            count+=1
            if f(a)*f(c)<0:
                b=c
            elif f(c)*f(b)<0:
                a=c
            elif f(c)==0:
                return [c,count]
                break
            c=(a+b)/2
        return [c,count]
    elif f(a)*f(b)==0:
        if f(a)==0:
            return a
        else:
            return b
#defining newton-raphson method
def NR_rt(F,f,x_start,x_tol):
    count=0
    x_prev=100000
    while abs(x_start-x_prev)>x_tol:
        count+=1
        x_prev=x_start
        x_start =x_start - F(x_start)/f(x_start)
    
    return [x_start, count]

#defining secant method
def sec_rt(f,a,b,x_tol):
    c=0
    count=0
    while abs(a-b)>=x_tol:
        count+=1
        c=b
        k=(f(b)- f(a))
        b = (f(b)*a - f(a)*b)/k
        a=c
        count+=1
    return [b,count]

def fn(x):
    return m.cos(x)-0.5

def dfn(x):
    return -m.sin(x)


tol=0.000000000001
print("The root found using NR method is {}, with {} number of iteration".format(NR_rt(fn,dfn,1.2,tol)[0],NR_rt(fn,dfn,1.2,tol)[1]))
print("The root found using secant method is {}, with {} number of iteration".format(sec_rt(fn,0.8,1.2,tol)[0],sec_rt(fn,0.8,1.2,tol)[1]))
print("The root found using bisection method is {}, with {} number of iteration".format(bisec_rt(fn,0.8,1.2,tol)[0],bisec_rt(fn,0.8,1.2,tol)[1]))        



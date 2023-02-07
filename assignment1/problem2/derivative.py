import numpy as np
import matplotlib.pyplot as plt

def dF2c(f,x,h):
    return (f(x+h)-f(x))/h

def dF3c(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)


#this part will generate different values of h
h_list=[1]

for i in range(0,20,2):
    h_list.append(h_list[i]/2)
    h_list.append(h_list[i]/10)

h_array=np.array(h_list)

#this part will show the value of derivative with different values of h
print("h                     dF2c                           dF3c")
x=2**0.5
dF2c_list=[]
dF3c_list=[]
for h in h_list:
    print(h,"                 ",dF2c(np.arctan,x,h),"                ",dF3c(np.arctan,x,h))
    dF2c_list.append(dF2c(np.arctan,x,h))
    dF3c_list.append(dF3c(np.arctan,x,h))

#error calculation
f2c_array=np.array(dF2c_list)
f3c_array=np.array(dF3c_list)

list_1_3=np.ones(len(f2c_array))/3

rel_err_2c= abs((f2c_array - list_1_3)*3)
rel_err_3c= abs((f3c_array - list_1_3)*3)

for i in range(len(h_list)):
    print(np.log(h_list[i]),"   ",np.log(rel_err_2c[i]),"     ",np.log(rel_err_3c[i]))

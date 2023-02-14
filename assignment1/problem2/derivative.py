import numpy as np
import matplotlib.pyplot as plt

#defining the two types of derivative
#2point derivative
def dF2c(f,x,h):
    return (f(x+h)-f(x))/h

#3point derivative
def dF3c(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)


#this part will generate different values of h
h_list=[1]

for i in range(0,40,2):
    h_list.append(h_list[i]/2)
    h_list.append(h_list[i]/10)

h_array=np.array(h_list)

x=np.sqrt(2) #taking the derivative at perticular x vlue

#part b
print('PART B: h vs single precision defivative and double precision derivative respectively:')
for h in h_list:
    print('{:.4e}       {:.4f}      {:.4f}'.format(h,dF2c(np.arctan,x,h),dF3c(np.arctan,x,h)))

print("\n\n\n\n")

#this part will show the value of derivative with different values of h
print('PART B: h vs single precision defivative and double precision derivative respectively:')
print("h                 dF2c         dF3c")
x=2**0.5
dF2c_list=[]
dF3c_list=[]
for h in h_list:
    print('{:.4e}       {:.4f}      {:.4f}'.format(h,dF2c(np.arctan,x,h),dF3c(np.arctan,x,h)))
    dF2c_list.append(dF2c(np.arctan,x,h))
    dF3c_list.append(dF3c(np.arctan,x,h))

print('\n\n\n\n')

#part c
#error calculation
f2c_array=np.array(dF2c_list)
f3c_array=np.array(dF3c_list)

list_1_3=np.ones(len(f2c_array))/3

rel_err_2c= abs((f2c_array - list_1_3)*3)
rel_err_3c= abs((f3c_array - list_1_3)*3)

print('PART C: h vs error of two method respectively')
for i in range(len(h_list)):
    print('{:.4e}       {:.4f}      {:.4f}'.format(np.log10(h_list[i]),np.log10(rel_err_2c[i]),np.log10(rel_err_3c[i])))


#plottings
#part b
#plotting the the derivative of different method with exact value (x axis in log scale)
plt.title('PART B: The value of derivative with different h')
plt.plot(h_array,dF2c_list,label='df2c')
plt.plot(h_array,dF3c_list,label='df3c')
plt.plot(h_array,list_1_3,label='Exact derivative')
plt.ylabel('Derivative')
plt.xlabel('h')
plt.xscale('log')
plt.legend()
plt.grid()
plt.show()

#part c
plt.title('PART C: log10(h) vs relative error')
plt.plot(np.log10(h_array),np.log10(rel_err_2c),label="rel err 2c")
plt.plot(np.log10(h_array),np.log10(rel_err_3c),label="rel err 3c")
plt.xlabel('log10(h)')
plt.ylabel('relative error')
plt.legend()
plt.grid()
plt.show()


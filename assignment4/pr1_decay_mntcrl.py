import numpy as np
import matplotlib.pyplot as plt

import numpy as np

def rad_decay(N0:int, lambda_, dt, Nruns):
    #initial conditions
    Nt= N0
    
    t_list_decay=[]
    Nt_list= []  #NbyM array
    t_list = []   #nbymarray
    for i in range(Nruns):
        count=0
        Nt_l=[]
        t_l= []
        while Nt>0 :
            p= lambda_ *dt*Nt #probability
            #generating random number b/w 0 & 1
            rand_x = np.random.rand()

            if rand_x <= p:
                Nt-=1   # one particle decay if rand_x satasfy the probability condition
            count+=1
            Nt_l.append(Nt)
            t_l.append(count*dt)
        Nt_list.append(Nt_l)
        t_list.append(t_l)
    
    return t_list, Nt_list

# Parameters
N0_values = [10, 100, 1000] # Initial number of nuclei
lambda_ = 0.3 # Decay constant
dt = 0.5* 10**(-3) # Time step
Nruns = 1000 # Number of runs/events

for N0 in N0_values:
    t_list, Nt_list = rad_decay(N0,lambda_,dt,Nruns)
    max_Nt= max(Nt_list,key=len)
    max_t= max(t_list,key=len)

    # Pad the sublists with zeros
    for sublist in Nt_list:
        sublist.extend([0] * (len(max_Nt) - len(sublist)))

    Nt_list= np.array(Nt_list)  #converting it to numpy array

    Nt_avg=np.zeros(len(max_Nt))
    for sl in Nt_list:
        Nt_avg+= sl
    
    Nt_avg =Nt_avg/Nruns

    plt.plot(max_t,Nt_avg)
plt.show()

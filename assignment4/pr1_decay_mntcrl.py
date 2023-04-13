import numpy as np
import matplotlib.pyplot as plt

import numpy as np

#function that do the simulation and give averege over Nruns
# N0 - number of particle at t = 0
# lambda_ - decay const
# dt - time steps
# Nruns - number of simulation
def rad_decay(N0:int, lambda_, dt, Nruns):
    #initial conditions
    Nt= N0
    
    t_list_decay=[]
    Nt_list= []  #NbyM array  #list of lists
    t_list = []   #nbymarray
    #running the test forNruns time
    for i in range(Nruns):
        count=0
        Nt_l=[]
        t_l= []
        #loop computing one decay
        while Nt>0 :
            p= lambda_ *dt*Nt #probability
            #generating random number b/w 0 & 1
            rand_x = np.random.rand()

            if rand_x <= p:
                Nt-=1   # one particle decay if rand_x satasfy the probability condition
            count+=1
            Nt_l.append(Nt)
            t_l.append(count*dt)
        Nt=N0
        Nt_list.append(Nt_l)
        t_list.append(t_l)
    
    return t_list, Nt_list   #lists of lists of time and N(t) for Nruns number of time

# Parameters
N0_values = [10, 100, 1000] # Initial number of nuclei
lambda_ = 0.3 # Decay constant
dt = 0.5* 10**(-3) # Time step
Nruns = 1000 # Number of runs/events


#running the different N0 value
for N0 in N0_values:
    t_list, Nt_list = rad_decay(N0,lambda_,dt,Nruns)
    #finding the list with maximum length(i.e the event which took the maximum time)
    max_Nt= max(Nt_list,key=len) 
    max_t= max(t_list,key=len)
    # Pad the sublists with zeros
    for sublist in Nt_list:
        sublist.extend([0] * (len(max_Nt) - len(sublist)))

    Nt_list= np.array(Nt_list)  #converting it to numpy array
    #compution the Nexact
    n_exact= N0* np.exp(-lambda_*np.array(max_t))
    #computing avg
    Nt_avg=np.zeros(len(max_Nt))
    for sl in Nt_list:
        Nt_avg+= sl  #summing up
    
    Nt_avg =Nt_avg/Nruns #taking < N(t) >
    log_nt_avg=np.log(Nt_avg)
    plt.plot(max_t,log_nt_avg, label='Decay of {} particles simulated'.format(N0))
    plt.plot(max_t,np.log(n_exact),ls=':', label='EXACT calculated result for {} particles'.format(N0))

    #calculating the derivative of the curve ln<N(t)> at t=0, O(h**2) accuracy
    df= -(log_nt_avg[0]-log_nt_avg[len(log_nt_avg)//5])/(dt*(len(log_nt_avg)//5))  # Second-order forward difference formula (accuracy O(h^2))
    print('the value of the derivative for N(0)={} is {}'.format(N0,df))

plt.xlabel('time (s)')
plt.ylabel('$ln<N(t)>$')
plt.title("Plot of time vs $lnN(t)$")
plt.legend()
plt.show()

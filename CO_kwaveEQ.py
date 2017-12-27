# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:28 2017

@author: JSitters
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#state Variables
model_segs=pd.DataFrame(data={'name':["BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6"],
                   'length':[9703, 6641, 5707, 9209, 5420, 10658], 
                   'bot_width':[1.22, 2.1699, 2.660, 2.830, 3.8599, 3.9199],
                   'chan_slope':[0.001760, 0.00171, 0.00063, 0.00163, 0.00028, 0.001],
                   'mannings_n':[0.06, 0.12, 0.075, 0.115, 0.065, 0.12]})
    
x=model_segs.length
n=model_segs.mannings_n
s=model_segs.chan_slope
w=model_segs.bot_width
beta=3/5
alpha=((n*w**(2/3))/s**0.5)**beta
delta_t=360
end_time=8514000
#setup


time2d=(np.arange(0,end_time,delta_t))
Qboun2d=(1-np.cos((time2d*np.pi)/450))*(750/np.pi)+200 #np.atleast_2d
Qx=np.ones((len(time2d),len(x)+2),dtype=float) 
Qout_CO=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)
Qout_Chap=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)
tt=[]


for i in range(x.count()):
    init_time=0
    for j in range(len(time2d)-1):
        this_time = Qout_CO[j][0]
        dt=this_time - init_time
        Qout_CO[j+1][i+2]=(delta_t/x[i]*Qout_CO[j+1][i+1] + alpha[i]*beta*Qout_CO[j][i+2]*((Qout_CO[j][i+1] + Qout_CO[j+1][i+1])/2)**(beta-1))/(delta_t/x[i] + alpha[i]*beta*((Qout_CO[j][i+2]+Qout_CO[j+1][i+1])/2)**(beta-1))
        Qout_Chap[j+1][i+2]=Qout_Chap[j][i+2]+((Qout_Chap[j][i+1]-Qout_Chap[j][i+2])*(dt))/(x[i]*alpha[i]*beta*Qout_Chap[j][i+2]**(beta-1))
        init_time = this_time
    
    tt.append(dt)

#%%     Plotting Code     #######   

plt.plot(Qout_CO[:,1], 'b', label='Boundary Flow')
plt.plot(Qout_CO[:,2], 'c', label='FirstCO')
plt.plot(Qout_Chap[:,2], 'g', label='FirstChap')
#plt.plot(Qout_CO[4], 'r', label='Third')
#plt.plot(Qout_CO[5], 'k', label='Fourth')
plt.ylabel('Flow')
plt.title('Colorado and Chapra Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,50,0,600])
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 11:39:28 2017

@author: JSitters
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

#state Variables
model_segs=pd.DataFrame(data={'name':["BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6"],
#                  'length':[9703, 6641, 5707, 9209, 5420, 10658], 
                                             'length':[900, 900,900, 900,900, 900], 
  
                   'bot_width':[1.22, 1.22, 1.22, 1.22, 1.22, 1.22],
                   'chan_slope':[0.001760, 0.001760, 0.001760,0.001760, 0.001760, 0.001760],
                   'mannings_n':[0.12,.12, 0.12, 0.12, 0.12, 0.12]})
                        #mannings_n':[0.06, 0.06, 0.06, 0.06, 0.06, 0.06]})
    
x=model_segs.length
n=model_segs.mannings_n
s=model_segs.chan_slope
w=model_segs.bot_width
beta=3/5
alpha=((n*w**(2/3))/s**0.5)**beta
delta_t=60
end_time=10140

#setup


time2d=(np.arange(0,end_time,delta_t))
Qboun2d=np.full((len(time2d)),5)
#(1-np.cos((time2d*np.pi)/450))*(750/np.pi)+2 #np.atleast_2d
Qx=np.ones((len(time2d),len(x)),dtype=float) 
Qout_CO=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)
Qout_Chap=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)
Qout_Chap[0,2:]=0.001*len(x)
Qout_CO[0,2:]=0.001*len(x)

for t in range(len(time2d)-1):
   
    for i in range(x.count()):
        
        Qout_CO[t+1][i+2]=((delta_t/x[i])*Qout_CO[t+1][i+1] + alpha[i]*beta*Qout_CO[t][i+2]*((Qout_CO[t][i+2] + Qout_CO[t+1][i+1])/2)**(beta-1))/(delta_t/x[i] + alpha[i]*beta*((Qout_CO[t][i+2]+Qout_CO[t+1][i+1])/2)**(beta-1))
        Qout_Chap[t+1][i+2]=Qout_Chap[t][i+2]+((Qout_Chap[t][i+1]-Qout_Chap[t][i+2])*(delta_t))/(x[i]*alpha[i]*beta*Qout_Chap[t][i+2]**(beta-1))


#%%    WASP input    
wasptime=range(0, 170,1)
os.getcwd()
os.chdir('C:/Users/jsitters/WASP/K_wave_Trial/')
WASP_Output = pd.read_table('WASP_200F.txt', sep="\t", header=0) #this is not right... needs work on segment location
WASP_df = pd.DataFrame(WASP_Output)  

#%%     Plotting Code     #######   
plt.figure(1)
#plt.plot(time2d/delta_t,Qout_CO[:,1], 'b', label='Boundary Flow')
plt.plot(wasptime,WASP_df.iloc[:,1], 'r', label ='WASP_seg1')
plt.plot(time2d/delta_t,Qout_CO[:,2], 'c', label='CO_seg1')
plt.plot(time2d/delta_t,Qout_Chap[:,2], 'g', label='Chap_seg1')
plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

plt.figure(2)
plt.plot(wasptime,WASP_df.iloc[:,2], 'r', label ='WASP_seg2')
plt.plot(time2d/delta_t,Qout_CO[:,3], 'c', label='CO_seg2')
plt.plot(time2d/delta_t,Qout_Chap[:,3], 'g', label='Chap_seg2')
plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

plt.figure(3)
plt.plot(wasptime,WASP_df.iloc[:,3], 'r', label ='WASP_seg3')
plt.plot(time2d/delta_t,Qout_CO[:,4], 'c', label='CO_seg3')

plt.plot(time2d/delta_t,Qout_Chap[:,4], 'g', label='Chap_seg3')
plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

plt.figure(4)
plt.plot(wasptime,WASP_df.iloc[:,4], 'r', label ='WASP_seg4')
plt.plot(time2d/delta_t,Qout_CO[:,5], 'c', label='CO_seg4')
plt.plot(time2d/delta_t,Qout_Chap[:,5], 'g', label='Chap_seg4')
plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

plt.figure(5)
plt.plot(wasptime,WASP_df.iloc[:,5], 'r', label ='WASP_seg5')
plt.plot(time2d/delta_t,Qout_CO[:,6], 'c', label='CO_seg5')
plt.plot(time2d/delta_t,Qout_Chap[:,6], 'g', label='Chap_seg5')
plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

plt.figure(6)
plt.plot(wasptime,WASP_df.iloc[:,6], 'r', label ='WASP_seg6')
plt.plot(time2d/delta_t,Qout_CO[:,7], 'c', label='CO_seg6')
plt.plot(time2d/delta_t,Qout_Chap[:,7], 'g', label='Chap_seg6')

plt.ylabel('Flow')
plt.xlabel('Time (min)')
plt.title('Kwave Comparison')
plt.legend(loc='lower right')
plt.axis([0,100,0,5.5])

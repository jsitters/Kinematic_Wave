# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:21:03 2017

@author: tSitters
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#state Variables
model_segs=pd.DataFrame(data={'name':["BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6"],
                   'length':[9703, 6641, 5707, 9209, 5420, 10658], 
                   'bot_width':[1.22, 2.1699, 2.660, 2.830, 3.8599, 3.9199],
                   'chan_slope':[0.001760, 0.00171, 0.00063, 0.00163, 0.00028, 0.001],
                   'mannings_n':[0.06, 0.12, 0.075, 0.115, 0.065, 0.12], 
                   'initial_depth':[1.73, 2.3, 2.41, 5.58, 2.8, 3.06]})
    
x=model_segs.length
n=model_segs.mannings_n
s=model_segs.chan_slope
w=model_segs.bot_width
delta_t=360
end_time=514000#8514000

##Channel Geometry 
beta=3/5
alpha=2

##setup

time2d=(np.arange(0,end_time,delta_t))
Qboun2d=(1-np.cos((time2d*np.pi)/450))*(750/np.pi)+200 
Qx=np.zeros((len(time2d),len(x)+2),dtype=float) 
Qx[0,:]=[0.001]*(len(x)+2) #initial flow conditions of 0.001 m3/s
Qout=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)

tt=[]


for t in range(len(time2d)-1):#switch these loops x on inside
    init_time=0
    this_time = Qout[t][0]
    dt=this_time - init_time
    for i in range(x.count()):
        #alpha=[((n/s**0.5)*(w+2*model_segs.initial_depth[i])**(2/3))**beta]
        Qout[t+1][i+2]=Qout[t][i+2]+((Qout[t][i+1]-Qout[t][i+2])*(dt))/(x[i]*alpha*beta*Qout[t][i+2]**(beta-1))
        init_time = this_time
        #area=alpha*(Qout[t+1][i+2])**beta
        #depth=area/w[i]
        #alpha=alpha.append(((n/s**0.5)*(w+2*depth)**(2/3))**beta)
        #velocity=Qout[t+1][i+2]/area
    tt.append(dt)
#%%
plt.plot(Qout[:100,0],Qout[:100,1], 'b', label='Boundary Flow')
plt.plot(Qout[:100,0],Qout[:100,2], 'c', label='First')
plt.plot(Qout[:100,0],Qout[:100,3], 'g', label='Second')
plt.plot(Qout[:100,0],Qout[:100,4], 'r', label='Third')
plt.plot(Qout[:100,0],Qout[:100,5], 'k', label='Fourth')
plt.ylabel('Flow')

#plt.axis([0,200,2000,2500])
#plt.legend() 
                
#%%
initial_depth=2 #meters initial depth
beta=3/5
z=2 #slope of 0.5

Perimeter=w+2*initial_depth #square
alpha1=((n/s**0.5)*(Perimeter)**(2/3))**beta
area=alpha1*Qout[t+1][i+2]**beta

#print(np.subtract(Qout[4][2],Qout[4][8])) #using numpy functions is faster than python functions


shape=['trapezoid', 'rectangle','triangle']
for i in shape:
    if shape=='trapezoid':
        Perimeter=w+2*initial_depth*(1+z**2)**0.5
    elif shape=='rectangle':
        Perimeter=2*initial_depth+w
    elif shape=='triangle':
        Perimeter=2*initial_depth*(1+z**2)
    alpha1=((n/s**0.5)*(Perimeter)**(2/3))**beta
print(alpha1)
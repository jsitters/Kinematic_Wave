# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:21:03 2017

@author: tSitters
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def perimeter_calc(shape, d,wide,z_slope):
    if shape =='tz':#for trapezoid
        perim=wide+2*d*(1+z_slope**2)**0.5    #this may not be calculating correctly
    elif shape =='r':#for rectangle
        perim=2*d+wide
    elif  shape =='t':#for triangle
        perim=2*d*(1+z_slope**2)**0.5
    return perim

#state Variables
model_segs=pd.DataFrame(data={'name':["BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6"],
                   'length':[9703, 6641, 5707, 9209, 5420, 10658], 
                   'bot_width':[1.22, 2.1699, 2.660, 2.830, 3.8599, 3.9199],
                   'chan_slope':[0.001760, 0.00171, 0.00063, 0.00163, 0.00028, 0.001],
                   'mannings_n':[0.06, 0.12, 0.075, 0.115, 0.065, 0.12], 
                   'initial_depth':[1.73, 2.3, 2.41, 5.58, 2.8, 3.06],
                   'chan_geom':['r', 'r', 'r', 'r', 'r','r']})
length=model_segs.length
n=model_segs.mannings_n
slope=model_segs.chan_slope
w=model_segs.bot_width
delta_t=60
end_time=600#8514000
model_segs.insert(loc= 7,column='Perimeter',value=0)#allow_duplicates=True


##setup dataframes

time2d=(np.arange(0,end_time,delta_t))
Qboun2d=time2d+200#(1-np.cos((time2d*np.pi)/450))*(750/np.pi)+200 #time2d
Qx=np.zeros((len(time2d),len(length)+1),dtype=float) 
Qx[0,:]=[0.001]*(len(length)+1)              #initial flow conditions of 0.001 m3/s
Qout=np.concatenate((Qboun2d[:,None],Qx), axis=1)

alpha_data=pd.DataFrame(data=np.zeros((len(time2d), length.count())), columns=list(model_segs.name))

area=pd.DataFrame(data=np.zeros((len(time2d), length.count())), columns=list(model_segs.name))
depth_data =pd.DataFrame(data=np.zeros((len(time2d), length.count())), columns=list(model_segs.name))
velocity=pd.DataFrame(data=np.zeros((len(time2d), length.count())), columns=list(model_segs.name))

##need to make the depth column in alpha_data be the initial depth. 
depth_data.loc[0,:]=list(model_segs['initial_depth'])


##Channel Geometry 
beta=3/5
z=2 #assuming trapezoid and triangle have a side slope of 1/2
#for j in range(len(model_segs['chan_geom'])):
#    if model_segs['chan_geom'][j] =='tz':#for trapezoid
#        model_segs['Perimeter'][j]=w[j]+1*depth_data[j]*(1+z**2)**0.5    #this may not be calculating correctly
#    elif model_segs['chan_geom'][j] =='r':#for rectangle
#        model_segs['Perimeter'][j]=2*depth_data[j]+w[j]
#    elif  model_segs['chan_geom'][j] =='t':#for triangle
#        model_segs['Perimeter'][j]=2*depth_data[j]*(1+z**2)**0.5
#print(model_segs.Perimeter)


##Time and Space Flow Loop####


for t in range(len(time2d)-1): 

    for i in range(length.count()):
        p_calc=perimeter_calc(model_segs['chan_geom'][i],depth_data.iloc[t,i], model_segs['bot_width'][i], z)
        
        #Perimeter=2*depth_data.iloc[t,i]+model_segs['bot_width'][i] #model_segs['Perimeter'][i]
        alpha_data.iloc[t,i]=((n[i]/slope[i]**0.5)*(p_calc**(2/3)))**beta
        Qout[t+1][i+1]=Qout[t][i+1]+((Qout[t][i]-Qout[t][i+1])*(delta_t))/(length[i]*alpha_data.iloc[t,i]*beta*Qout[t][i+1]**(beta-1))

        area.iloc[t,i]=alpha_data.iloc[t,i]*(Qout[t+1][i+1])**beta
        depth_data.iloc[t+1,i]=area.iloc[t,i]/w[i]
        #alpha_data.iloc[t,i]=alpha
        velocity.iloc[t,i]=Qout[t+1][i+1]/area.iloc[t,i]
        
print(alpha_data.head(n=6))
print(depth_data.head(n=6))
print(area.head(n=6))        
    
#print(alpha_data[0:15]) # shows zeros for the last time step

#%%
##Plotting Code##

plt.plot(time2d,Qout[:100,0], 'b', label='Boundary Flow')
plt.plot(time2d,Qout[:100,1], 'c', label='First')
plt.plot(time2d,Qout[:100,2], 'g', label='Second')
plt.plot(time2d,Qout[:100,3], 'r', label='Third')
plt.plot(time2d,Qout[:100,4], 'k', label='Fourth')
plt.ylabel('Flow')

#plt.axis([0,200,2000,2500])
plt.legend() 
                


#print(np.subtract(Qout[4][2],Qout[4][8])) #using numpy functions is faster than python functions





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
                   'initial_depth':[1.73, 2.3, 2.41, 5.58, 2.8, 3.06],
                   'chan_geom':['r', 'r', 'tz', 'tz', 't','t']})
x=model_segs.length
n=model_segs.mannings_n
s=model_segs.chan_slope
w=model_segs.bot_width
delta_t=60
end_time=86400#8514000
model_segs.insert(loc= 7,column='Perimeter',value=0)#allow_duplicates=True

#model_segs['alpha'] = ((model_segs.mannings_n*model_segs.bot_width**(2/3))/model_segs.chan_slope**0.5)**beta

##setup dataframes

time2d=(np.arange(0,end_time,delta_t))
Qboun2d=(1-np.cos((time2d*np.pi)/450))*(750/np.pi)+200 
Qx=np.zeros((len(time2d),len(x)+2),dtype=float) 
Qx[0,:]=[0.001]*(len(x)+2)              #initial flow conditions of 0.001 m3/s
Qout=np.concatenate((time2d[:,None],Qboun2d[:,None],Qx), axis=1)
alpha_data=pd.DataFrame(np.zeros(((len(time2d)*x.count()), x.count()-1)))
alpha_data.columns =['alpha', 'area', 'depth', 'velocity', 'width']

##need to make the depth column in alpha_data be the initial depth. 
alpha_data.loc[0:6,'depth']=model_segs['initial_depth']
##Channel Geometry 
beta=3/5
z=2 #assuming trapezoid and triangle have a side slope of 1/2
for j in range(len(model_segs['chan_geom'])):
    if model_segs['chan_geom'][j] =='tz':#for trapezoid
        model_segs['Perimeter'][j]=w[j]+2*alpha_data['depth'][j]*(1+z**2)**0.5    #this may not be calculating correctly
    elif model_segs['chan_geom'][j] =='r':#for rectangle
        model_segs['Perimeter'][j]=2*alpha_data['depth'][j]+w[j]
    elif  model_segs['chan_geom'][j] =='t':#for triangle
        model_segs['Perimeter'][j]=2*alpha_data['depth'][j]*(1+z**2)**0.5
#print(model_segs.Perimeter)


##Time and Space Flow Loop####
f = 0
for t in range(len(time2d)-1):
    
    for i in range(x.count()):
        #model_segs['chan_geom'][i]
        alpha=((n[i]/s[i]**0.5)*(model_segs['Perimeter'][i])**(2/3))**beta
        Qout[t+1][i+2]=Qout[t][i+2]+((Qout[t][i+1]-Qout[t][i+2])*(delta_t))/(x[i]*alpha*beta*Qout[t][i+2]**(beta-1))

        alpha_data['area'][i+f]=alpha*(Qout[t+1][i+2])**beta
        alpha_data['depth'][i+f]=alpha_data['area'][i+f]/w[i]
        alpha_data['alpha'][i+f]=alpha
        alpha_data['velocity'][i+f]=Qout[t+1][i+2]/alpha_data['area'][i+f]
        alpha_data['width'][i+f]=w[i]
        
    f = f + x.count()
print(alpha_data[8600:]) # shows zeros for the last time step
alpha_data.shape
#%%
##Plotting Code##

plt.plot(Qout[:100,0],Qout[:100,1], 'b', label='Boundary Flow')
plt.plot(Qout[:100,0],Qout[:100,2], 'c', label='First')
plt.plot(Qout[:100,0],Qout[:100,3], 'g', label='Second')
plt.plot(Qout[:100,0],Qout[:100,4], 'r', label='Third')
plt.plot(Qout[:100,0],Qout[:100,5], 'k', label='Fourth')
plt.ylabel('Flow')

#plt.axis([0,200,2000,2500])
plt.legend() 
                


#print(np.subtract(Qout[4][2],Qout[4][8])) #using numpy functions is faster than python functions





# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 11:30:42 2018

@author: JSitters
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def perimeter_calc(shape, d, wide, z_slope):
    if shape =='tz':#for trapezoid
        perim=wide +2 *d*(1 + z_slope**2)**0.5    
    elif shape =='r':#for rectangle
        perim=2*d + wide
    elif  shape =='t':#for triangle
        perim=2 *d*(1 + z_slope**2)**0.5
    return perim

#state Variables
beta=3/5
delta_t=3600/2
endtime=8514000/6
time_s=np.arange(0, endtime, delta_t)
#flow_in = [200]*5 + [10]*5 + [50]*5 + [100]*75 + [100]*4640 #Boundary flow
     
model_segs = pd.DataFrame(data={'name':["BC_1", "BC_2", "BC_3", "BC_4", "BC_5", "BC_6"],
                   'length':[9703, 6641, 5707, 9209, 5420, 10658], 
                   'bot_width':[29.62, 29.62, 29.62, 29.62, 29.62, 29.62],
                   'chan_slope':[0.001760, 0.00171, 0.00063, 0.00163, 0.00028, 0.001],
                   'mannings_n':[0.06, 0.12, 0.075, 0.115, 0.065, 0.12],
                   'init_depth':[1.73, 2.3, 2.41, 2.58, 2.8, 3.06],
                   'init_vol':[20397, 33083, 36551, 67245, 58605, 127845],          ##never used
                   'ChanGeom':["r", "r", "r", "r", "r", "r"],
                   'z_slope':[2, 2, 2, 2, 2, 2],   
                   'depth_exp':[0.5343, 0.5343, 0.5343, 0.5343, 0.5343, 0.5343],
                   'depth_intercept':[-0.8234, -0.8234, -0.8234, -0.8234, -0.8234, -0.8234]})
    

##output setup
nsegments=len(model_segs)
names=list(model_segs['name'])
names.insert(0,"t")
names.insert(1,"flow_in")
Qout=pd.DataFrame(np.zeros((len(time_s)+1,nsegments+2)), columns=names)
Qout.iloc[0,2:]=list([200]*(nsegments))              #initial flow conditions of 0.001 m3/s

## Hydrologic Parameters data setup
alpha_data=pd.DataFrame(data=np.zeros((len(time_s), nsegments)), columns=list(model_segs.name))
HP_vol=pd.DataFrame(data=np.zeros((len(time_s), nsegments+2)), columns=names)
HP_area=pd.DataFrame(data=np.zeros((len(time_s), nsegments+2)), columns=names)
HP_velocity=pd.DataFrame(data=np.zeros((len(time_s), nsegments+2)), columns=names)
HP_depth =pd.DataFrame(data=np.zeros((len(time_s)+1, nsegments+2)), columns=names)
HP_depth.loc[0,2:]=list(model_segs['init_depth'])       #initital depth plugged into the depth datafram

##add first two columns to dataframes
Qout['t'] = np.concatenate((time_s,[time_s.max() + delta_t]))
HP_depth['t']=np.concatenate((time_s,[time_s.max() + delta_t]))
HP_vol['t']=time_s
HP_velocity['t']=time_s    
HP_area['t']=time_s    
Qout['flow_in'] = [200]*10 + [10]*10 + [50]*10 + [100]*(len(time_s)-29)      ##User Defined Boundary flow condition
HP_depth['flow_in'] = [200]*10 + [10]*10 + [50]*10 + [100]*(len(time_s)-29)
HP_velocity['flow_in'] = [200]*10 + [10]*10 + [50]*10 + [100]*(len(time_s)-30)  
HP_area['flow_in']=[200]*10 + [10]*10 + [50]*10 + [100]*(len(time_s)-30)
HP_vol['flow_in']=[200]*10 + [10]*10 + [50]*10 + [100]*(len(time_s)-30)
#pd.options.mode.chained_assignment = None  # default='warn'

##Loop through segments to calculate outflow and hydro parameters for each time step
b=0
for t in range(len(time_s)): 
    for i in range(nsegments):
        p_calc=perimeter_calc(model_segs['ChanGeom'][i],HP_depth.iloc[t,i+2], model_segs['bot_width'][i], model_segs['z_slope'][i])
        alpha_data.iloc[t,i]=((model_segs['mannings_n'][i]/model_segs['chan_slope'][i]**0.5)*(p_calc**(2/3)))**beta

        Q_hat = Qout.iloc[t, i+2]
        #Qout.iloc[t+1, i+2]= Qout.iloc[t, i+2]+((Qout.iloc[t, i+1]-Qout.iloc[t, i+2])*(delta_t))/(model_segs['length'][i]*alpha_data.iloc[t,i]*beta*Qout.iloc[t, i+2]**(beta-1))
        Qout.iloc[t+1, i+2] = Qout.iloc[t, i+2] + ((Qout.iloc[t+1, i+1] - Qout.iloc[t, i+2])*(delta_t))/(model_segs['length'][i]*alpha_data.iloc[t, i]*beta*Q_hat**(beta-1))
        while Q_hat - Qout.iloc[t+1, i+2] > 1e-3:  # implicit soln for Qout.iloc[t+1, i+2]
            Q_hat = Qout.iloc[t+1, i+2]
            Qout.iloc[t+1, i+2] = Qout.iloc[t, i+2] + ((Qout.iloc[t+1, i+1] - Qout.iloc[t, i+2])*(delta_t))/(model_segs['length'][i]*alpha_data.iloc[t, i]*beta*Q_hat**(beta-1))
            b = b + 1
            print(b)
        #else Q_hat - Qout.iloc[t+1, i+2] <= 1e-3:
            #Qout.iloc[t+1, i+2] = Qout.iloc[t, i+2] + ((Qout.iloc[t+1, i+1] - Qout.iloc[t, i+2])*(delta_t))/(model_segs['length'][i]*alpha_data.iloc[t, i]*beta*Q_hat**(beta-1))
        HP_area.iloc[t,i+2]=alpha_data.iloc[t,i]*(Qout.iloc[t+1 , i+2])**beta
        HP_depth.iloc[t+1,i+2]=HP_area.iloc[t,i+2]/model_segs['bot_width'][i]
        HP_velocity.iloc[t,i+2]=Qout.iloc[t+1 , i+2]/HP_area.iloc[t,i+2]
        HP_vol.iloc[t,i+2]=model_segs['length'][i]*HP_area.iloc[t,i+2]
     
        ## may need to remove first initial conditions row from HP_depth and Qout  
HP_Output=pd.concat([HP_depth, HP_area, HP_velocity, HP_vol], keys=['depth', 'area', 'velocity', 'volume'])

os.chdir('C:/Users/jsitters/WASP/K_wave_Trial/Time Step/')

#np.savetxt("WASP_equ.csv", Qout, delimiter=',')  

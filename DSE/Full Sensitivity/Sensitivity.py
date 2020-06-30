# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 14:27:19 2020

@author: Burhan
"""

#import Range_sense as ranger
import numpy as np
from matplotlib import pyplot as plt
import math

g = 9.81
M = 0.78
a = 295.0423
t_loiter = 1800
R_divert = 500000
PLW = 20000

LD_cruise = 28.11
V_cruise = 230.1
LD_loiter = 29.22
TSFC_cruise = 1.14e-5
TSFC_loiter = 0.8*TSFC_cruise

MTOW = 70469
OEW  = 41272
MF   = 9188

a = 0.431
b = 10903   #slop and constant 


Mass_Ratio = MTOW/OEW
LD_ratio = LD_cruise/LD_loiter

Fuel_cap = 20000


def F_loiter(LD_l = LD_loiter, TSFC_l = TSFC_loiter):
    
    
    F_l = math.exp(-t_loiter/(LD_l/(g*TSFC_l)))
    return F_l
    
def F_divert(LD_cruise, TSFC_crusie):
    
    F_d = math.exp(-R_divert/(LD_cruise*V_cruise/(g*TSFC_crusie)))    
    return F_d

def cruise_range(LD_c, TSFC_c, w_f):
    
    c_range = (V_cruise/(TSFC_c*g))*LD_c*math.log(w_f)
    
    return (c_range)


def find_prev_mass(M_current, F_current):
    
    '''
    Function to find fuel mass left at the end of various post-cruise phases
    based on the mass at the end of the curent psot-cruise phase
    
    Inputs:
        M_current = mass at the end of the current post-cruise phase
        F_current = fuel mass fraction at the end of the current post cruise phase
        
    Outputs:
        M_prev = mass at the end of the previous post-cruise phase
    '''
    
    M_prev = (M_current/F_current) + (1/F_current - 1)*(OEW + PLW)
    return M_prev

def find_next_mass(M_current, F_next):
    
    '''
    Function to find the mass at the end of the next pre-cruise phase based
    on the mass at the end of the current post cruise phase
    
    Inputs:
        M_current = mass at the end of the current pre-crusie phase
        F_next    = fuel fraction of the next pre-cruise phase
        
    Outputs:
        M_next    = mass at the end of the next cruise phase
    '''
    
    M_next = F_next*M_current - (1 - F_next)*(OEW + PLW)   
    return M_next


def find_weight_frac(OEW, LD_c, TSFC_c, LD_l, TSFC_l):
    F_pre_cruise = np.array([0.99, 0.99, 0.995, 0.98])
    F_post_cruise = np.array([F_divert(LD_c, TSFC_c), F_loiter(LD_l, TSFC_l), 0.99, 0.992])
    
    M_array = np.zeros(10)
    M_array[0] = MF
    
    for i in range (len(F_post_cruise)):
        
        M_array[-i-2] = find_prev_mass(M_array[-i-1], F_post_cruise[-i-1])
        
    for j in range (len(F_pre_cruise)):
        
        M_array[j+1] = find_next_mass(M_array[j], F_pre_cruise[j])
        
    # f_precruise = 1 - M_array[4]/M_array[0]
    # f_postcruise = 1 - M_array[7]/M_array[0]
    # f_reserve = (M_array[5] - M_array[7])/M_array[0]
    # f_cruise = (M_array[4] - M_array[5])/M_array[0]
    
    wf = (OEW + PLW + M_array[4])/(OEW + PLW + M_array[5])
    
    return wf

    

OEW_arr = np.arange(0.9*OEW, 1.1*OEW, 10)
MF_arr  = np.zeros(len(OEW_arr))
MTOW_arr = np.zeros(len(OEW_arr))
LD_arr = np.linspace(0.9*LD_cruise, 1.1*LD_cruise, len(OEW_arr))
TSFC_arr = np.linspace(1.1*TSFC_cruise, 0.9*TSFC_cruise, len(OEW_arr))
range_arr = np.zeros(len(OEW_arr))
efficiency_arr = np.zeros(len(OEW_arr))

OEW_parr = (OEW_arr - OEW)/OEW*100


for k in range (len(OEW_arr)):

    OEW = OEW_arr[k]
    #LD_cruise = LD_arr[k]
    #TSFC_cruise = TSFC_arr[k]
    
    TSFC_loiter = 0.8*TSFC_cruise
    LD_loiter = LD_cruise/LD_ratio
    
    
    weight_frac = find_weight_frac(OEW, LD_cruise, TSFC_cruise, LD_loiter, TSFC_loiter)
    
    range_C = cruise_range(LD_cruise, TSFC_cruise, weight_frac)
    
    while range_C/1000 + 400 > 3995.5 or range_C/1000 + 400 < 3995:
    
        while range_C/1000 + 400 < 3995:
            #print ('Hello')
            MF+=0.2
            MTOW+=0.2
            #OEW = a*MTOW + b
            MTOW = PLW + OEW + MF
            
            weight_frac = find_weight_frac(OEW, LD_cruise, TSFC_cruise, LD_loiter, TSFC_loiter)
            range_C = cruise_range(LD_cruise, TSFC_cruise, weight_frac)
            
        while range_C/1000 + 400> 3995.5:
            #print ('Bye')
            MF-=0.2
            MTOW-=0.2
            #OEW = a*MTOW + b
            MTOW = PLW + OEW + MF
            
            weight_frac = find_weight_frac(OEW, LD_cruise, TSFC_cruise, LD_loiter, TSFC_loiter)
            range_C = cruise_range(LD_cruise, TSFC_cruise, weight_frac)
        
    range_arr[k] = range_C/1000 + 400
    
    efficiency_arr[k] = MF/3995/194
    
    print (k)
        
    MF_arr[k] = MF
    MTOW_arr[k] = MTOW
    
    total_range = 400 + range_C/1000    
    range_arr[k] = total_range
    
    
range_parr = (range_arr - 3995)/3995*100
MTOW_parr = (MTOW_arr - 70469)/70469*100
MF_parr = (MF_arr - 9188)/9188*100
LD_parr = (LD_arr - 28.11)/28.11*100
TSFC_parr = (TSFC_arr - 1.14e-5)/1.14e-5*100
efficiency_parr = (efficiency_arr - 0.01185)/0.01185*100

plt.rcParams.update({'font.size': 18})
fig = plt.figure(figsize = (10,10))
plt.plot(OEW_parr, efficiency_arr, label = "CNA fuel burn")
plt.grid(True)
plt.hlines(0.01431, -10, 10, colors = 'g', label = "A320neo fuel burn")
plt.hlines(0.01431*0.9, -10, 10, colors = 'r', label = "A320neo fuel burn reduced by 10%")
plt.xlabel("Change in L/D cruise [%]")
plt.ylabel(r"Fuel burn $\left[\frac{kg}{pax \cdot km}\right]$")
plt.legend(loc = "upper left")
plt.ylim([0.9*min(efficiency_arr), 1.1*0.01431])
plt.show() 
    
    
    
    





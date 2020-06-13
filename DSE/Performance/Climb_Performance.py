#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:58:13 2020

@author: youssef
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sp
# INPUTS
MTOW = 73593.64093*9.81        # maximum take off weight [N]
S    = 150.7210057            # wing surface area [m^2]


g    = 9.81                   # gravitational acceleration [m/s^2]
y    = 1.4
R    = 287 
T0   = 288.15
TSFC = 1.12*10**-5    # kg/Ns

Tto  = 189873.8014             # thrust at take off [N]

BPR  = 15                     # bypass ratio
rho0 = 1.225                  # air density surface level [kg/m^3]
AR   = 17
e    = 0.7
k1   = 0.046881 
CD0  = 0.02005
CLmax_clean = 1.7536
# parameters 


# functions

def ISA_trop(h):
	""" This function computes the atmospheric properties 
		within the troposphere
	Input:
		h = altitude [m]
	Output:
		T = Temperature [K]
		p = Pressure [Pa]
		rho = Densituy [kg/m^3]
		a = spped of sound [m/s]
	"""
	T = 288.15 - 0.0065*h;
	p = 101325*(T/288.15)**(-g/(-0.0065*287));
	rho = 1.225*(T/288.15)**(-g/(-0.0065*287) - 1);
	a = np.sqrt(1.4*287*T);
	return T, p, rho, a;

def Net_Thrust(V,h):
    
    """

    Parameters
    ----------
    M : mach number [ - ]

    Returns
    -------
    Tnet = net thrust [N]

    """
    M = V/(ISA_trop(h)[3])
    n = 0.7
    sigma = (ISA_trop(h)[2])/rho0
    if M <= 0.4:
        K1 = 1.0
        K2 = 0.0
        K3 = -0.595
        K4 = -0.03
    if M > 0.4:
        K1 = 0.89
        K2 = -0.014
        K3 = -0.300
        K4 = 0.005
    

    Tnet = Tto*((K1 + K2*BPR + (K3 +K4*BPR)*M)*(sigma**n))
    return Tnet

# STEADY FLIGHT CONDITIONS, PLOT POWER VS VELOCITY, sea level condition

   # m/s

def ROC_st_f(h):
    V  = list(range(50,300,1))
    Pa  = []
    Pr  = []
    Vr  = []
    for i in range(len(V)):
        
        CL = 2*MTOW/((ISA_trop(h)[2])*S*(V[i]**2))
        if CL >= CLmax_clean:
            continue
        
        Ma = V[i]/(ISA_trop(h)[3])
        if Ma >= 0.78:
            continue
        CD = CD0 + k1*(CL**2)
        D  = 1/2*(ISA_trop(h)[2])*(V[i]**2)*CD*S
        T  = Net_Thrust(V[i],h)
        
        Pa.append(T*V[i])
        Pr.append(D*V[i])
        Vr.append(V[i])
    Pa = np.array(Pa)
    Pr = np.array(Pr)
    ROC_st = (Pa-Pr)/MTOW*196.85
    plt.plot(Vr,ROC_st)
    return ROC_st

<<<<<<< HEAD
altitude = np.arange(0,35000,20000)
=======






def ROC_unst_f(h):
    V    = list(range(50,300,1))
    Pa   = []
    Pr   = []
    Vr   = []
    DVDH = []
    for i in range(len(V)):
        
        CL = 2*MTOW/((ISA_trop(h)[2])*S*(V[i]**2))
        # check CLmax_clean is respected
        if CL >= CLmax_clean:
            continue
        Ma = V[i]/(ISA_trop(h)[3])
        #check maximum mach number is respected
        if Ma >= 0.78:
            continue
        CD = CD0 + k1*(CL**2)
        D  = 1/2*(ISA_trop(h)[2])*(V[i]**2)*CD*S
        T  = Net_Thrust(V[i],h)
        # correction for unsteady flight
        Veas = V[i]*np.sqrt((ISA_trop(h)[2])/rho0)
        drhodH = (y/T0)*((g/(2*R*y))+1/2)*((1+y*h/T0)**(g/(2*R*y)-1/2))
        dVdH = Veas*drhodH
        #append values to list
        Pa.append(T*V[i])
        Pr.append(D*V[i])
        Vr.append(V[i])
        DVDH.append(dVdH)
    Pa       = np.array(Pa)
    Pr       = np.array(Pr)
    Vr       = np.array(Vr)
    DVDH     = np.array(DVDH)
    ROC_st   = (Pa-Pr)/MTOW
    ROC_unst = ROC_st/(1+Vr*DVDH/g)*196.85
    plt.plot(Vr,ROC_unst)
    return ROC_unst

'''
altitude = np.arange(0,40000,5000)
>>>>>>> 92af3e55e454a052bcc51b7911b02db3a2b9b71c

for j in range(len(altitude)):
    ROC_st_f(0.3048*altitude[j])



    
plt.grid(True)
plt.ylabel('ROC [ft/min]')
plt.xlabel('Velocity [m/s]')    
plt.show()
'''



# Thrust function with interpolation

Altitude = np.arange(0,13000,1000)
Mach     = np.array([0,0.25,0.5,0.78])
Values   = np.array([[96722.5216	,   67428.05362,47797.17295	,34848.43582],
                     [92730.66919,  66953.61247, 49363.17223,34875.76088],
                     [83327.33083,	60272.56754,	44580.51898,	31811.88102],
                     [75528.07162,	54984.87745,	41105.11362,	30166.20097],
                     [68192.30809,	49935.16789,	37672.74605,	28178.95313],
                     [61325.46673,	45143.27947,	34329.63267,	26076.81658],
                     [54934.35083,	40632.33221,	31115.11061,	23940.46769],
                     [52728.40771,	38788.01903,	29339.68686,	21994.39766],
                     [47022.7086,	34789.833,	26546.78776,	20260.27162],
                     [41496.00717,	30791.69771,	23600.46648,	18162.6575],
                     [36575.19359,	27245.49901,	21000.56856,	16329.92886],
                     [33957.59999,	25411.71708,	19713.45868,	15500.42165],
                     [27471.22224,	20554.11247,	15941.42122,	12529.85862]])

Thrustf = sp.interp2d(Mach,Altitude,Values,kind='cubic')
step = 5000 #ft
altitude = np.arange(0,40000,step)*0.3048

def ROC_st_true(h,MTOW):
    
    
    #h = array containing all altitudes
    

    Machs = np.arange(0.2,0.78,0.01)
    
    Thrust = Thrustf(Machs,h)   # rows=altitudes; columns= mach numbers
    for i in range(len(Thrust)):
        
        height = h[i]   
        V = []
        for j in range(len(Machs)):
            Vm = ISA_trop(height)[3]*Machs[j]
            V.append(Vm)
        
       
        Pa=[]
        Pr=[]
        Vr=[]
        Tr=[]
        for k in range(len(V)):
        
            CL = 2*MTOW/((ISA_trop(height)[2])*S*(V[k]**2))
            
            if CL >= CLmax_clean:
                continue
            CD = CD0 + k1*(CL**2)
            
            D  = 1/2*(ISA_trop(height)[2])*(V[k]**2)*CD*S
            
            T  = 2*Thrust[i,k]
            
            Pa.append(T*V[k])
            Pr.append(D*V[k])
            Vr.append(V[k])
            Tr.append(T)
        Pa = np.array(Pa)
        Pr = np.array(Pr)
        Tr = np.array(Tr)
        ROC_st = (Pa-Pr)/MTOW*196.85
        '''
        if np.amax(ROC_st) < 0:
            plt.plot(Vr,ROC_st)
            continue
        '''
        
        FW = TSFC*Tr[np.where(ROC_st == np.amax(ROC_st))[0]][0]*step*0.3048/(np.amax(ROC_st)/196.85)
        MTOW -= FW
        print(MTOW,FW)
        
        plt.plot(Vr,ROC_st)
    

ROC_st_true(altitude,MTOW)
plt.grid()
plt.xlabel("Velocity [m/s]")
plt.ylabel("ROC steady [ft/min]")
plt.show()

















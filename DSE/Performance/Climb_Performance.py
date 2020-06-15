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
MTOW = (79584.74877 - 1500)*9.81        # maximum take off weight [N]
MTOWi = MTOW
S    = 150.865002            # wing surface area [m^2]
OEW  = 45818.75425*9.81
PW   = 20000*9.81
mcruise = 76000*9.81
g    = 9.81                   # gravitational acceleration [m/s^2]
y    = 1.4
R    = 287 
T0   = 288.15
TSFc = 1.203380285E-05  # kg/Ns

Tto  = 211576.8505             # thrust at take off [N]

BPR  = 15                     # bypass ratio
rho0 = 1.225                  # air density surface level [kg/m^3]
AR   = 17
e    = 0.7
k1   = 0.044763
k2   = 1/np.pi/AR/e
CD0  = 0.017423
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
# Thrust function with interpolation

Altitude = np.arange(0,13000,1000)
Mach     = np.array([0,0.25,0.5,0.78])
Values   = 1*np.array([[130983.30,	93336.08,67645.75,	50789.30],
                     [122348.90,	87785.60,	64908.05	,44054.64],
                     [112232.70,	79860.39,	57123.37,	43527.24],
                     [101170.00,	72339.87,	55721.66,	41455.83],
                     [91730.72,	66149.63	,48671.44,	35633.27],
                     [82740.84,	60087.01,	44755.67	,33719.74],
                     [74338.93,	54325.67,	40861.95,	31365.01],
                     [72199.00,	52709.15,	39372.54,	29609.87],
                     [64465.65,	47372.79	, 35782.72,	27582.68],
                     [57326.40,	42373.03,	32291.17,	25290.77],
                     [50687.43,	37660.26	,29101.37,	23136.07],
                     [44974.74,	33673.71	,26141.89,	21070.06],
                     [38392.90,	28740.28,	22306.38,	18085.59]])

Thrustf = sp.interp2d(Mach,Altitude,Values,kind='linear')
# TSFC function with interpolation
Altitude = np.arange(0,16000,4000)
Mach     = np.arange(0,1,0.5)
Values   = (10**-6)*np.array([[5.1996, 10.07987],
                     [5.53191581,	9.747087486],
                     [6.117945716,	10.29400869],
                     [6.655300979,	10.95611142]])

TSFCf    = sp.interp2d(Mach,Altitude,Values,kind='linear')

step = 1000*0.3048 #ft
altitude = np.arange(0*0.3048,40000*0.3048,step)

def ROC_st_true(h,MTOW):
    
    
    #h = array containing all altitudes
    
    HD = []
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
        DVDH = []
        for k in range(len(V)):
        
            CL = 2*MTOW/((ISA_trop(height)[2])*S*(V[k]**2))
            
            if CL >= CLmax_clean:
                continue
            CD = CD0 + k1*(CL**2)
            
            D  = 1/2*(ISA_trop(height)[2])*(V[k]**2)*CD*S
            
            T  = 2*Thrust[i,k]
            # correction for unsteady flight
            Veas = V[k]*np.sqrt((ISA_trop(height)[2])/rho0)
            drhodH = (y/T0)*((g/(2*R*y))+1/2)*((1+y*height/T0)**(g/(2*R*y)-1/2))
            dVdH = Veas*drhodH
            Pa.append(T*V[k])
            Pr.append(D*V[k])
            Vr.append(V[k])
            Tr.append(T)
            DVDH.append(dVdH)
        Pa = np.array(Pa)
        Pr = np.array(Pr)
        Tr = np.array(Tr)
        Vr = np.array(Vr)
        DVDH = np.array(DVDH)
        
        ROC_st = (Pa-Pr)/MTOW*196.85
        #ROC_st = np.divide(ROC_st,(1+1/g*np.multiply(Vr,DVDH)))
        
        if np.amax(ROC_st) < 0:
            plt.plot(Vr,ROC_st)
            continue
        
        
        
        #horizontal distance
        Vd   = Vr[np.where(ROC_st == np.amax(ROC_st))[0]][0]
        #print(Vd)
        gamma = np.arcsin(np.amax(ROC_st)/196.85/Vd)
        hd = (step/np.tan(gamma))*0.001
        
        
        HD.append(hd)
    
        # fuel weight reduction
        
        
        Macr          = Vr/ISA_trop(height)[3]    # mach numbers for 
        TSFC          = TSFCf(Macr,h)[i]
        TSFCrocmax    = TSFC[np.where(ROC_st == np.amax(ROC_st))[0]][0]
        Trocmax       = Tr[np.where(ROC_st == np.amax(ROC_st))[0]][0]
        time_to_climb = step/(np.amax(ROC_st)/196.85)
        FW  = 9.81*TSFCrocmax*Trocmax*time_to_climb
        
        MTOW = MTOW - FW
        #print(Tr[np.where(ROC_st == np.amax(ROC_st))[0]][0] ,np.amax(ROC_st) )
        '''
        if height == 30000*0.3048:
            print(np.amax(ROC_st))
        '''
        
        plt.plot(Vr,ROC_st)
    print("horizontal distance  =",sum(HD)-HD[-1],"km")
    print("fuel consumed        =",(MTOWi - MTOW)/9.81, "kg")
    

ROC_st_true(altitude,MTOW)

plt.grid()
plt.xlabel("Velocity [m/s]")
plt.ylabel("ROC steady [ft/min]")
plt.show()

















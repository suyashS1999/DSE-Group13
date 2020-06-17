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
MTOW = (70459 - 1402 - 345)*9.81        # maximum take off weight [N]
MTOWi = MTOW
S    = 150.865002            # wing surface area [m^2]
OEW  = 41271.3770*9.81
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
k2   = 1/np.pi/AR/e

k1   = 0.039648
CD0  = 0.009312


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

Altitude = np.arange(0,15000,3000)
Mach     = np.array([0,0.25,0.5,0.75])

Values   = 1*np.array([[83846.38395,	58099.14016,	40336.57846,	30877.3181],
                        [64975.86428,	46995.09957,	34718.0064,	25653.11629],
                        [47418.72899,	34910.376,	26530.64582,	20692.9504],
                        [37554.91053,	27962.75161,	21446.6353,	16886.62497],
                        [26253.48281,	19691.85709,	15223.2655,	12116.85906]])

Thrustf = sp.interp2d(Mach,Altitude,Values,kind='linear')
# TSFC function with interpolation
Altitude = np.arange(0,15000,3000)
Mach     = np.array([0,0.25,0.5,0.75])


Values   = (10**-6)*np.array([[5.078422255,	7.199381363,	10.1043708,	13.13544401],
                                [5.464774141,	7.444791219,	9.627267933,	12.01378245],
                                [5.896717801,	7.909076521,	10.01071562,	11.98768685],
                                [5.676529015,	7.542071193,	9.514016207,	11.40673502],
                                [5.674835113,	7.49283021,	9.409177574,	11.22866788]])


TSFCf    = sp.interp2d(Mach,Altitude,Values,kind='linear')

step = 100*0.3048 #ft
altitude = np.arange(0,40000*0.3048,step)

MRC_arr = []
VMRC_arr = []
Steps_arr = []
HD_Arr = []

def ROC_st_true(h,MTOW):
    
    
    #h = array containing all altitudes

    hdarr = 0

    checked = 0
    
    HD = []
    Machs = np.arange(0.2,0.75,0.01)
    
    Thrust = Thrustf(Machs,h)   # rows=altitudes; columns= mach numbers

    print (len(Thrust))
    
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
        #DVDH = []
        for k in range(len(V)):
        
            CL = 2*MTOW/((ISA_trop(height)[2])*S*(V[k]**2))
            
            if CL >= CLmax_clean:
                continue
            CD = CD0 + k1*(CL**2)
            
            D  = 1/2*(ISA_trop(height)[2])*(V[k]**2)*CD*S
            
            T  = 2*Thrust[i,k]
            # correction for unsteady flight
            #Veas = V[k]*np.sqrt((ISA_trop(height)[2])/rho0)
            #drhodH = (y/T0)*((g/(2*R*y))+1/2)*((1+y*height/T0)**(g/(2*R*y)-1/2))
            #dVdH = Veas*drhodH
            Pa.append(T*V[k])
            Pr.append(D*V[k])
            Vr.append(V[k])
            Tr.append(T)
            #DVDH.append(dVdH)
        Pa = np.array(Pa)
        Pr = np.array(Pr)
        Tr = np.array(Tr)
        Vr = np.array(Vr)
        #DVDH = np.array(DVDH)
        
        ROC_st = (Pa-Pr)/MTOW*196.85 #ft/min
        #ROC_st = np.divide(ROC_st,(1+1/g*np.multiply(Vr,DVDH)))
        
        if np.amax(ROC_st) < 500 and checked == 0:
            plt.plot(Vr,ROC_st, 'o')
            print("The ceilling is - ", h[i]/0.3048, "ft")
            checked = 1
            continue
        
        maxroc = np.argmax(ROC_st)

        MAXROC = ROC_st[maxroc]

        MRC_arr.append(MAXROC)
        
        
        #horizontal distance
        Vd   = Vr[maxroc]
        VMRC_arr.append(Vd)
        #print(Vd)
        gamma = np.arcsin(np.amax(ROC_st)/196.85/Vd)

        #V_hor = Vr[maxroc]*np.cos(gamma)


        hd = (step/np.tan(gamma))*0.001

        hdarr+=hd
        
        
        HD.append(hd)

        HD_Arr.append(hdarr)
    
        # fuel weight reduction
        
        
        Macr          = Vr/ISA_trop(height)[3]    # mach numbers for
        TSFC          = TSFCf(Macr,h)[i]
        TSFCrocmax    = TSFC[maxroc]
        Trocmax       = Tr[maxroc]
        time_to_climb = step/(np.amax(ROC_st)/196.85)
        FW  = 9.81*TSFCrocmax*Trocmax*time_to_climb
        
        MTOW = MTOW - FW
        #print(Tr[np.where(ROC_st == np.amax(ROC_st))[0]][0] ,np.amax(ROC_st) )
        '''
        if height == 30000*0.3048:
            print(np.amax(ROC_st))
        '''
        
        #plt.plot(Vr,ROC_st)
    #plt.plot(VMRC_arr, MRC_arr)
    
    HD_Arr.append(HD_Arr[-1])
    print (len(HD_Arr), len(altitude))
    plt.plot(HD_Arr, altitude)
    print("horizontal distance  =",sum(HD)-HD[-1],"km")
    print("fuel consumed        =",(MTOWi - MTOW)/9.81, "kg")
    

ROC_st_true(altitude,MTOW)

plt.grid()
plt.xlabel("Distance [km]")
plt.ylabel("Altitude [m]")
#plt.xlabel("Velocity [m/s]")
#plt.ylabel("ROC steady [ft/min]")
plt.show()

















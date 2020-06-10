#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:58:13 2020

@author: youssef
"""
import numpy as np
import matplotlib.pyplot as plt
# INPUTS
MTOW = 74616.9829*9.81        # maximum take off weight [N]
S    = 139.1091985            # wing surface area [m^2]
g    = 9.81                   # gravitational acceleration [m/s^2]

Tto  = 211545.862             # thrust at take off [N]

BPR  = 15                     # bypass ratio
rho0 = 1.225                  # air density surface level [kg/m^3]
AR   = 17
e    = 0.7
k    = 1/(np.pi*AR*e)
CD0  = 0.01632
CLmax_clean = 1.112
# parameters 

G = 0.06*BPR + 0.64   # High BPR engines  https://books.google.it/books?id=ZQu5DwAAQBAJ&pg=PA205&lpg=PA205&dq=empirical+formula+for+thrust+to+altitude&source=bl&ots=crHjDLskeW&sig=ACfU3U3LTfDSQEHYZHqOFy_HkBlkpAsF3g&hl=it&sa=X&ved=2ahUKEwjJyuu16vbpAhXCCuwKHZg-CYAQ6AEwDXoECAoQAQ#v=onepage&q=empirical%20formula%20for%20thrust%20to%20altitude&f=false
k1 = 0.377*(1+BPR)/np.sqrt((1+0.82*BPR)*G)
k2 = 0.23 +0.19*np.sqrt(BPR)


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

def Net_Thrust2(V,h):
    """

    Parameters
    ----------
    M : mach number [ - ]

    Returns
    -------
    Tnet = net thrust [N]

    """
    M = V/(ISA_trop(h)[3])
    Tnet = Tto*(1- k1*M +k2*M**2)
    return Tnet

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
    sigma = rho0/(ISA_trop(h)[2])
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
    

    Tnet = Tto*(K1 + K2*BPR + (K3 +K4*BPR)*M*(sigma**n))
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
        CD = CD0 + k*(CL**2)
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

altitude = np.arange(0,40000,5000)

for j in range(len(altitude)):

    ROC_st_f(0.3048*altitude[j])

plt.grid(True)
plt.title('Steady Rate of Climb')
plt.ylabel('ROC [ft/min]')
plt.xlabel('Velocity [m/s]')    
plt.show()



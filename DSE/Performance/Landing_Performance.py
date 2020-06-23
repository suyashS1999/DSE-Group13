#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:15:46 2020

@author: youssef
"""
import numpy as np
hscr       = 50*0.3048    # screen height [m]
G          = 9.81
mu_free    = 0.03  # friction coefficient free rolling (0.03-0.05)
mu_brake   =  0.3  # friction coefficient braking (0.3-0.5)
CLmaxland  = 2.4
MTOW       = 79584.74877*9.81        # maximum take off weight [N]
MLW        = 0.98*MTOW
rho        = 1.225 
g          = 9.81
S          = 150.865002
CD0clean   = 0.01204                # zero lift drag coefficient
# estimate drag contribution due to flap deflection

Rf = 0.25   # flap chord/wing chord
Swet_to_Sref = 0.65  # ratio between wetted area and wing surface area
Sflap_to_Sref =Rf*Swet_to_Sref
defl = 50 #degrees

delta1 = 179.32*Rf**4 - 111.6*Rf**3 + 28.929*Rf**2 + 2.3705*Rf - 0.0089
delta2 = (-3.9877*10**-12)*defl**6 + (1.1685*10**-9)*defl**5 + (-1.2846*10**-7)*defl**4 + (6.1742*10**-6)*defl**3 +(-9.89444*10**-5)*defl**2 + (6.8324*10**-4)*defl**-4 + (-3.892*10**-4)

deltaCDflap = delta1*delta2*Sflap_to_Sref


# ground effect

k    =  0.039648         # 1/(pi*AR*eto) from simulation
b    = 50.61874255      # wing span [m]
df   = 4.2              # fuselage diameter [m]
hlg  = 1.9              # clearance between surface and bottom fuselage due to landing gear [m]
e    = 0.75              # oswald efficiency factor at takeoff [-]
AR   = 17               # Aspect Ratio [-]

#calculations
Vsland   = np.sqrt((2/rho)*(MLW/S)*(1/CLmaxland))






# APPROACH DISTANCE
Vapp   = 1.3*Vsland
theta_app = 3 #degrees

n   = 1/2*rho*(Vapp**2)*0.9*CLmaxland*S/MLW #1.521
R   = Vapp**2/(g*(n-1))
hf  = R*(1-np.cos(np.radians(theta_app)))

Sa  = (hscr-hf)/np.tan(np.radians(theta_app))



# FLARE DISTANCE
Vfr = Vapp

Sf  = R*np.sin(np.radians(theta_app))
# FREE-ROLL DISTANCE
Vtd = 1.1*Vsland

Sfr = 3*Vtd


# BRAKING DISTANCE
Vbr    = 1.1*Vsland
Vg     = Vbr/np.sqrt(2)
CLbr  = CLmaxland/(1.1**2)  
CDiOGE = k*(CLbr**2)
h      = df + hlg
ge     = 1 - (1-1.32*h/b)/(1.05+7.4*h/b)
CDiIGE = ge*CDiOGE
CDbr  = CD0clean + CDiIGE + deltaCDflap
Lg     = 1/2*rho*(Vg**2)*CLbr*S
Dg     = 1/2*rho*(Vg**2)*CDbr*S
Tg     = 0            # thrust at Vbr/sqrt(2) assumed to be zero
mugnd  = 0.3         # ground friction coefficient 0.03-0.05 brakes off, otherwise 0.3-0.5 with braking

Sg     = -(Vbr**2)*MLW/(2*g*(Tg-Dg-mu_brake*(MLW-Lg)))


Sland= Sg + Sfr + Sf + Sa

print("Sa = ",Sa," m")
print("Sf = ",Sf," m")
print("Sfr = ",Sfr," m")
print("Sg = ",Sg," m")
print("Sland = ",Sland," m")

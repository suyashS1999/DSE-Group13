#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:04:38 2020

@author: youssef
"""
import numpy as np
#                BALANCED FIELD LENGTH BFL (Torenbeek)
# definition: runway distance required to reach V1 from stop and then hard brake until stop
# inputs

MTOW = 79584.74877*9.81     # maximum take off weight [N]
S    = 150.865002           # wing surface area
rho  = 1.225                # air density [kg/m^3]
g    = 9.81                 # gravitational acceleration [m/s^2]
CLto = 2.1                  # Clmax at take off [-]
hto  = 11                   # obstacle height [m], 35 ft for jetliners
Tto  = 211576.8505          # thrust at take off [N]
y2m  = 0.024                # some weird parameter for 2 engines
BPR  = 15                   # bypass ratio
CD0takeoff = 0.019482       # drag coefficient at takeoff ASSUMED
Sto  = 200                  # inertia distance
rho0 = 1.225                # air density surface level [kg/m^3]
CL0  = 0.182587             # lift coefficient at zero angle of attack
k    = 0.044698             # 1/(pi*AR*eto) from simulation
b   = 50.61874255           # wing span [m]
df  = 4.37                  # fuselage diameter [m]


#balanced field length

# calculations


Vs   = np.sqrt((2/rho)*(MTOW/S)*(1/CLto))
VMCA = 1.2*Vs
V2   = 1.1*VMCA
CL2  = 0.694*CLto
Tavg = 0.75*Tto*(5+BPR)/(4+BPR)
mu   = 0.01*CLto + 0.02
D2   = 1/2*rho*(V2**2)*(CD0takeoff + k*(CLto**2))*S
TOEI = Tto/2
y2   = np.arcsin((TOEI-D2)/MTOW)
dy2  = y2 - y2m
sigma = rho/rho0

BFL = 0.863/(1+2.3*dy2)*((MTOW/S)/(rho*g*CL2)+hto)*(2.7 + 1/(Tavg/MTOW-mu))+(Sto/np.sqrt(sigma))


#                                   GROUND ROLL
#inputs


hlg = 2              # clearance between surface and bottom fuselage due to landing gear [m]


#calculations

Vlof   = 1.1*Vs
CLlof  = 2*MTOW/(rho*(Vlof**2)*S)  

# ground effect
  
CDiOGE = k*(CLlof**2)
h      = df + hlg
ge     = 1 - (1-1.32*h/b)/(1.05+7.4*h/b)
CDiIGE = ge*CDiOGE

# flap deflection effect

Rf = 0.25   # flap chord/wing chord
Swet_to_Sref = 0.65  # ratio between wetted area and wing surface area
Sflap_to_Sref =Rf*Swet_to_Sref
defl = 10 #degrees

delta1 = 179.32*Rf**4 - 111.6*Rf**3 + 28.929*Rf**2 + 2.3705*Rf - 0.0089
delta2 = (-3.9877*10**-12)*defl**6 + (1.1685*10**-9)*defl**5 + (-1.2846*10**-7)*defl**4 + (6.1742*10**-6)*defl**3 +(-9.89444*10**-5)*defl**2 + (6.8324*10**-4)*defl**-4 + (-3.892*10**-4)

deltaCDflap = delta1*delta2*Sflap_to_Sref

# merge all together
Vg     = Vlof/np.sqrt(2)
Lg     = 1/2*rho*(Vg**2)*CL0*S
CDg    = CD0takeoff + CDiIGE
Dg     = 1/2*rho*(Vg**2)*CDg*S
Tg     = 105013.7572*2   # thrust at VLOF/sqrt(2)
mugnd  = 0.05         # ground friction coefficient 0.03-0.05 brakes off, otherwise 0.3-0.5 with braking


Sg     = (Vlof**2)*MTOW/(2*g*(Tg-Dg-mugnd*(MTOW-Lg)))

print("Sg  = ", Sg, " m")

# rotation


Sr  = 3*Vlof     # large aircraft take between 2 to 5 seconds to rotate
print("Sr  = ", Sr, " m")
# transition

#1 average vertical acceleration in terms of weight factor

n    = 1.1903

Vtr  = (Vlof +V2)/2
nu   = (Vtr/Vs)**2
R    = (Vtr**2)/(g*(nu-1))
CLtr = 2*MTOW/(rho*(Vtr**2)*S)
CLv2 = CLto/(1.2**2)
CDv2 = CD0takeoff + k*CLv2**2
CDtr = 0.0247 + (CLtr**2)*k
Dtr  = MTOW*CDtr/CLtr
Dv2  = 1/2*rho*(V2**2)*CDv2*S
Tc   = 94049.13479*2    # thrust at V2  Tto*(1-2*V2/np.sqrt(1.4*287*288.15)*(1+BPR)/(3+2*BPR))
thetaclimb = np.arcsin((Tc-Dv2)/MTOW)   # radians
htr        = R*(1-np.cos(thetaclimb))
Str        = np.sqrt((R**2)-(R-htr)**2)

#print("Str = ",Str," m")

Sobst = np.sqrt((R**2)-(R-hto)**2)

print("Sobst = ",Sobst," m")

# overall takeoff distance

Sto = Sg + Sr + Sobst 

print("Sto = ",Sto," m")


# second climb with one engine inoperative

gammaOEI = np.arcsin((Tc/2-Dv2)/MTOW)

print(gammaOEI*180/np.pi)



















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

MTOW = 74616.9829*9.81        # maximum take off weight [N]
S    = 139.1091985            # wing surface area
rho  = 1.225                  # air density [kg/m^3]
g    = 9.81                   # gravitational acceleration [m/s^2]
CLto = 2.0                    # Clmax at take off [-]
hto  = 10.7                   # obstacle height [m], 35 ft for jetliners
Tto  = 211545.862             # thrust at take off [N]
y2m  = 0.024                  # some weird parameter for 2 engines
BPR  = 15                     # bypass ratio
CDto = 0.05                   # drag coefficient at takeoff ASSUMED
Sto  = 200                    # inertia distance
rho0 = 1.225                  # air density surface level [kg/m^3]

# calculations

Vs   = np.sqrt((2/rho)*(MTOW/S)*(1/CLto))
V2   = 1.2*Vs
CL2  = 0.694*CLto
Tavg = 0.75*Tto*(5+BPR)/(4+BPR)
mu   = 0.01*CLto + 0.02
D2   = 1/2*rho*(V2**2)*CDto*S
TOEI = Tto/2
y2   = np.arcsin((TOEI-D2)/MTOW)
dy2  = y2 - y2m
sigma = rho/rho0



BFL = 0.863/(1+2.3*dy2)*((MTOW/S)/(rho*g*CL2)+hto)*(2.7 + 1/(Tavg/MTOW-mu))+(Sto/np.sqrt(sigma))








#                                   GROUND ROLL
#inputs

b   = 48.62978895      # wing span [m]
df  = 4.2              # fuselage diameter [m]
hlg = 2.1              # clearance between surface and bottom fuselage due to landing gear [m]
eto = 0.75             # oswald efficiency factor at takeoff [-]
AR  = 17               # Aspect Ratio [-]


#calculations
CDiOGE = (CLto**2)/(np.pi*AR*eto)
Vlof   = 1.1*Vs
h      = df + hlg
ge     = 1 - (1-1.32*h/b)/(1.05+7.4*h/b)
CDiIGE = ge*CDiOGE
Llof   = 1/2*rho*((Vlof/np.sqrt(2))**2)*CLto*S
Dlof   = Llof*CDto/CLto
Tlof   = Tto
mugnd  = 0.05    # ground friction coefficient 0.03-0.05 brakes off, otherwise 0.3-0.5 with braking


Sg     = (Vlof**2)*MTOW/(2*g*(Tlof-Dlof-mugnd*(MTOW-Llof)))

print("Sg  = ", Sg, " m")

# rotation


Sr  = 3*Vlof     # large aircraft take between 2 to 5 seconds to rotate
print("Sr  = ", Sr, " m")
# transition

#1 average vertical acceleration in terms of weight factor

n    = 1.1903
Vtr  = 1.15*Vs
R    = (Vtr**2)/(g*(n-1))
V2   = 1.2*Vs
CLtr = 2*MTOW/(rho*(Vtr**2)*S)
CLv2 = CLto/(1.2**2)
CDv2 = 0.08399
CDtr = 0.0247 + (CLtr**2)/(np.pi*AR*eto)
Dtr  = MTOW*CDtr/CLtr
Dv2  = 1/2*rho*(V2**2)*CDv2*S
Tc   =  #Tto*(1-2*V2/np.sqrt(1.4*287*288.15)*(1+BPR)/(3+2*BPR))
thetaclimb = np.arcsin((Tc-Dtr)/MTOW)   # radians
htr        = R*(1-np.cos(thetaclimb))
Str        = R*np.sin(thetaclimb)

print("Str = ",Str," m")


# climb

Sc = (hto-htr)/np.tan(thetaclimb)


print("Sc = ",Sc," m")


# overall takeoff distance

Sto = Sg + Sr + Str + Sc

print("Sto = ",Sto," m")





















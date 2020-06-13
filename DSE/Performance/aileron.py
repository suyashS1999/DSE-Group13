#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:27:07 2020

@author: youssef
"""
import numpy as np

S         = 150.865002   # wing area [m2]
MTOW      = 79584.74877
rho0      = 1.225
CL_max_to = 2.1
Vs        = np.sqrt((2/rho0)*(MTOW*9.81/S)*(1/CL_max_to))
VMCA      = 1.2*Vs
V2               = 1.1*VMCA # V2 A320-200 [m/s]
b                = 18           # wing span inboard [m]
bi_b2            = 0.6          # inboard span, range 0.6-0.67   [-] ref
bo_b2            = 0.92         # outboard span, range 0.92-0.92 [-] ref
delta_a_up       = 25           # max upward deflection,   range 25-30 [deg] ref
delta_a_down     = 20           # max downward deflection, range 10-25 [deg] ref or 75% of delta_a_up
ca_c             = 0.4          # ratio chords, range 0.24-0.3

# inboard airfoil data

Clalpha = 6.446      # lift curve slope [1/rad]
Cd0     = 0.005       # zero lift coefficient NEEDS TO BE CHANGED
cr      = 3.973      # root chord [m]
cm      = 2.326      # mid chord  [m]
Sin     = 115.1      # surfcace area inboard [m^2]
tpin    = 0.59       # taper ratio inboard   [-]

# chord equation c(y)
a = 1/(b/(cm-cr))
r = cr
def c(y):
    c = a*y + r
    return c

# differential ailerons: determine max deflection for minimising adverse yaw effect

delta_a = np.radians(0.5*(delta_a_up + delta_a_down))    # radians


#                           CALCULATIONS

# b1 and b2

# b1 = bi_b2*b
# b2 = bo_b2*b

b1 = 16
b2 = 17

# integral ratio

I = (1/3*a*(b2**3-b1**3)+1/2*r*(b2**2-b1**2))/(1/4*a*b**4 + 1/3*r*b**3)

# Roll requirement

P = np.radians(60)/11    # 60 degrees roll in not more than 11 seconds [deg/s]

# actual roll

pi = 0.6      # aileron effectiveness (graph, depends on chord ratio)

Pr = Clalpha*pi/(Clalpha+Cd0)*I*delta_a*V2

print("P required =",P ) 
print("P actual   =",Pr)

if Pr >= P:
    print("Requirement satisfied? Yes ")
else:
    print("Requirement satisfied? No  ")







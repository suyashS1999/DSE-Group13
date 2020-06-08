# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:18:26 2020

@author: steph
"""

import matplotlib.pyplot as plt
from InterpolationandIntegration import M, V, y_m, Moment_dist
import math as ma

#print(M)
#print(V)
#print(y_m)

L = -V[-1]       #Total lift force in newtons
#print(L) 

ML = M[-1]
#print(ML)
#%%------------Generate simple lift distribution---------------------------
b_half = 24.3
l = L/b_half
#print(l)


#%%------------Equilibrium equations---------------------------
Fe = 2305.990439*9.81
e = 10
#print(Fe)
#M.append(-Fe*e)
#y_m.append(e)

d = 49.22/2*0.55

h = 18
F_H = 0.0981
#M.append(-F_H*h)
#y_m.append(h)

D = 4.2
gamma = ma.atan(d/D)

F_brz = ((ML-e*Fe-h*F_H)/d)
#M.append(-F_brz*d)
#y_m.append(d)

Ay = -F_brz/ma.tan(gamma)
Az = Fe - F_H + F_brz - L

#print(F_brz)

#%%-----------Bending axial stress calculation
E = 69000
I = 148500000*10**(-4)*10**(-4)*10**(-4)
y = 356.09999999999997/2*10**(-3)

#M_min
M_max = 341662.565  #Nm
x_M_max = 17.02 #m
sigma = M_max*y/I
print(sigma*10**(-6))

#%%-----------Bending maximum shear stress
t = 2*10**(-3)              #mm
h = 356.09999999999997*10**(-3)     #mm
w = 1068.3*10**(-3)                 #mm

y = (h/2*t*h/4*2 + t*(w - 2*t)*(h/2 - t/2))/(h/2*t*2 + t*(w-2*t))  #m

A = (h/2*t*2 + t*(w - 2*t))

Q = A*y
V_max = 130837.069

tau = V_max*Q/(I*t)
print(tau*10**(-6))

sigma_vonmises = ma.sqrt((sigma**2 + 6*tau**2)/2)
print(sigma_vonmises*10**(-6))
sigma_yield = 324*10**6

if sigma_vonmises < sigma_yield:
    print("okay")
    
else: 
    print("not okay")
#print(sigma*10**(-6))
#print(sigma_yield*10**(-6))

#%%-------------------Wingbox mass estimation---------
b_half_inboard = 18
Fuel_vol = b_half_inboard*(w-2*t)*(h-2*t) #m^3
print("Fuel Space per wing", Fuel_vol)

rho_al = 2780
Mat_vol = h*w*b_half - (h-2*t)*(w-2*t)*b_half #m^3
print("Wingbox material volume =", Mat_vol)
Wingmass = Mat_vol*rho_al
print("Wingbox mass per side =", Wingmass)

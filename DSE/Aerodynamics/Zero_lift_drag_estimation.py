# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:02:17 2020

@author: sanjay
"""



import numpy as np
from Aerodynamics_module2 import *


# Wing Parameters

S_wet_ratio = 2.14         #- 
b_total = 49.22         #m
t_c_avg  = 0.15         # Aver. t/c
x_c_m    = 0.1          # location of max thickness
C_r        = 5          #m
C_t     = 2             #m
t_c_stream = 0.1        #-
CL_cruise  = 0.4744     #-
M_crit_airfoil = 0.592  #-   

sweep_LE = np.radians(40.63) 
sweep_max_thickness = sweep_x(x_c_m,sweep_LE,C_r, C_t, b_total)
M_dd_cruise = M_dd(t_c_stream, sweep_LE, CL_cruise)

rho_cruise = 0.3636     # kg/m^3
V_cruise   = 230.13     # m/s
mu_cruise = 1.46e-5     #SI units

laminar_flow = 0.1      #%
turb_flow = 1-laminar_flow  #%

k_wing = 0.634e-5        # Assuming smooth paint on wing

M_cruise = 0.78 

Re_cruise_max = Re_transonic(rho_cruise, V_cruise, C_r,mu_cruise, k_wing, M_cruise)


Cf_total = laminar_flow*Cf_laminar(Re_cruise_max) + turb_flow*Cf_turb(M_cruise, Re_cruise_max)

FF_wing = Form_factor_wing(t_c_avg, x_c_m, M_cruise, sweep_max_thickness)

M_crit_wing = M_crit_airfoil/np.cos(sweep_LE)

if M_crit_wing>M_cruise:
    delta_cd = 0
    print("Gooood")

elif M_crit_wing<=M_cruise and M_cruise<=M_dd_cruise:
    delta_cd = 0.002*((1+(2.5*((M_dd_cruise-M_cruise)/0.05)))**-1)
    print("Okeii")
elif M_dd_cruise<M_cruise:
    delta_cd = 0.002*((1+(((M_dd_cruise-M_cruise)/0.05)))**25)
    print("Baaaaaaddd")
    
C_d_misc = 5        #% of total CD0


CD_0_wing = (1+(C_d_misc/100))*((S_wet_ratio*Cf_total*FF_wing) + delta_cd)

print("The total Drag (CD0) : {} ".format(CD_0_wing))





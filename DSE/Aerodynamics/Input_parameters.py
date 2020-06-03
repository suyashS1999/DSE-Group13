# -*- coding: utf-8 -*-
"""

"""

import numpy as np




AR      = 17            #-
b_total = 49.22         #m
t_c_avg  = 0.15         # Aver. t/c
x_c_m    = 0.1          # location of max thickness
C_r     = 5             #m
C_t     = 2             #m
t_c_stream = 0.1        #-
CL_cruise  = 0.4744     #-
M_crit_airfoil = 0.592  #-  
S_wet_ratio = 2.14         #- 

sweep_LE = np.radians(40.63) 


# e  = e(AR,sweep_LE)
e = 0.7

rho_cruise = 0.3636     # kg/m^3
V_cruise   = 230.13     # m/s
mu_cruise = 1.46e-5     #SI units

laminar_flow = 0.1      #%
turb_flow = 1-laminar_flow  #%

k_wing = 0.634e-5        # Assuming smooth paint on wing

M_cruise = 0.78 
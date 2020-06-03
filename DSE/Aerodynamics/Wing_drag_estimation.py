# -*- coding: utf-8 -*-
"""

"""



import numpy as np
import Aero_tools as tl
from Input_parameters import *




"""Lift Calculations"""









"""Drag Calculation"""


CD_misc = 5        #% of total CD0

CD_0 = (1+(C_d_misc/100))*tl.CD0_wing(t_c_avg,x_c_m,sweep_LE,C_r,C_t,b_total,rho_cruise,V_cruise,mu_cruise,k_wing,M_cruise)

delta_CD = tl.delta_wave_drag(M_cruise, t_c_stream, sweep_LE, CL_cruise, M_crit_airfoil)

CD_induced = (CL_cruise**2)/(np.pi*AR*e)

Total_CD = CD_0 + CD_induced + delta_CD

print("The total Drag (CD) : {} ".format(Total_CD))





# -*- coding: utf-8 -*-
"""

"""



import numpy as np
import Aero_tools as tl
from Input_parameters import *


sweep_half_c = planf.sweep_x(0.5, sweep_LE, C_r_m, C_t_m, span)



"""Lift Calculations"""

#INBOARD


CL_alp_inb = tl.CL_alpha_DATCOM(AR_inb, M_cruise, a0_inb, sweep_half_c)

#OUTBOARD

CL_alp_outb = tl.CL_alpha_DATCOM(AR_outb, M_cruise, a0_outb, sweep_half_c)





"""Drag Calculation"""


CD_misc = 5        #% of total CD0

CD_0_cruise_inb = (1+(CD_misc/100))*tl.CD0_wing(t_c_avg_inb,x_c_m_inb,sweep_LE,C_r_inb,C_t_inb,span_inb,rho_cruise,V_cruise,mu_cruise,k_wing,M_cruise)

delta_CD = tl.delta_wave_drag(M_cruise, t_c_stream, sweep_LE, CL_cruise, M_crit_airfoil)

CD_induced = (CL_cruise**2)/(np.pi*AR*e)

Total_CD = CD_0 + CD_induced + delta_CD

print("The total Drag (CD) : {} ".format(Total_CD))





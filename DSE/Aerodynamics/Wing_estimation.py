# -*- coding: utf-8 -*-
"""

"""



import numpy as np
import Aero_tools as tl
from Input_parameters import *
import time as time

t_s = time.time()

sweep_half_c = tl.sweep_x(0.5, sweep_LE, C_r_m, C_t_m, span)



"""Lift Calculations"""

#INBOARD


CL_alp_inb = tl.CL_alpha_DATCOM(AR_inb, M_cruise, a0_inb, sweep_half_c)

#OUTBOARD

CL_alp_outb = tl.CL_alpha_DATCOM(AR_outb, M_cruise, a0_outb, sweep_half_c)





"""Drag Calculation"""


CD_misc = 5        #% of total CD0

#INBOARD

CD_0_cruise_inb = (1+(CD_misc/100))*tl.CD0_wing(t_c_avg_inb,x_c_m_inb,sweep_LE,C_r_m,C_t_m,Cr_inb,span,rho_cruise,V_cruise,mu_cruise,k_wing,M_cruise)

delta_CD_inb = tl.delta_wave_drag(M_cruise, t_c_stream_inb, sweep_LE, CL_cruise, M_crit_airfoil_inb)

CD_induced_inb = (CL_cruise**2)/(np.pi*AR_inb*e)

Total_CD_inb = CD_0_cruise + delta_CD_inb + CD_induced_inb

print("The total Drag (CD) INBOARD : {} ".format(Total_CD))


#OUTBOARD

CD_0_cruise_outb = (1+(CD_misc/100))*tl.CD0_wing(t_c_avg_outb,x_c_m_outb,sweep_LE,C_r_m,C_t_m,Cr_outb,span,rho_cruise,V_cruise,mu_cruise,k_wing,M_cruise)

delta_CD_outb = tl.delta_wave_drag(M_cruise, t_c_stream_outb, sweep_LE, CL_cruise, M_crit_airfoil_outb)

CD_induced_outb = (CL_cruise**2)/(np.pi*AR_outb*e)

Total_CD_outb = CD_0_cruise_outb + delta_CD_outb + CD_induced_outb


t_final = time.time()



"""Print Results"""


print("Inboard Section Results ..... \n \n")
print("CL_alpha  =  \n")
print("Total_CD  =  \n")
print("CL_max    =   \n ")
print("Stall angle =   \n \n \n ")

print("Outboard Section Results ..... \n \n")
print("CL_alpha  =  \n")
print("Total_CD  =  \n")
print("CL_max    =   \n")
print("Stall angle =   \n \n \n ")

print("Time Taken:", t_final-t_s, "seconds")





# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import Aero_tools as tl

"""Input Parameters"""

M_crit_wing = 0.78
M_cruise = 0.78

# Horizontal Tail

a0_h = 6.195564526   #1/rad
S_h = 35.03211872    #m^2
CL_h_design = 0.5
M_crit_a_h = tl.bisection(tl.M_cr_calc_h,0.2,0.99,100)
# print(M_crit_a_h)
A_h =   1           # Design Choice
taper_h = 0.85      # Design Choice
sweep_le_h = np.radians(42)     # deg

# Vertical Tail

a0_v = 6.195564526     #1/rad
S_v =  18.4     # m^2
CL_v_design = 0.1613575552564497
M_crit_a_v = tl.bisection(tl.M_cr_calc_v,0.2,0.99,100)
# print(M_crit_a_v)
A_v =   1         # Design Choice
taper_v = 0.85     # Design Choice
sweep_le_v = np.radians(42)    # deg


"""Computations"""

# Horizontal Tail

b_h = np.sqrt(A_h*S_h)
Cr_h = (2*S_h)/(b_h*(1+taper_h))
Ct_h = taper_h*Cr_h
sweep_half_h = tl.sweep_x(0.5,sweep_le_h,Cr_h,Ct_h,b_h)
CL_h_alpha = tl.CL_alpha_DATCOM(A_h,M_cruise,a0_h,sweep_half_h)
M_crit_h = M_crit_a_h/np.cos(sweep_le_h)

if M_crit_h>=(M_crit_wing+0.05):
	val_h = "Good"
else:
	val_h = "Bad"

# Vertical Tail

b_v = np.sqrt(A_v*S_v)
Cr_v = (2*S_v)/(b_v*(1+taper_v))
Ct_v = taper_v*Cr_v
sweep_half_v = tl.sweep_x(0.5,sweep_le_v,Cr_v,Ct_v,b_v)
CL_v_alpha = tl.CL_alpha_DATCOM(A_v,M_cruise,a0_v,sweep_half_v)
M_crit_v = M_crit_a_v/np.cos(sweep_le_v)

if M_crit_v>=(M_crit_wing+0.05):
	val_v = "Good"
else:
	val_v = "Bad"

"""Print Results"""

print("Horizontal Tail \n \n")
print("Span =", b_h)
print("CL_alpha =", CL_h_alpha)
print("Cr =", Cr_h)
print("Ct =",Ct_h)
print("M_crit =",M_crit_h, val_h,"\n \n")


print("Vertical Tail \n \n")
print("Span =", b_v)
print("CL_alpha =", CL_v_alpha)
print("Cr =", Cr_v)
print("Ct =",Ct_v)
print("M_crit =",M_crit_v, val_v, "\n \n")

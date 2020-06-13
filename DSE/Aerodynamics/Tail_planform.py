# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import Aero_tools as tl

"""Input Parameters"""

M_crit_wing = 0.78
M_cruise = 0.78

# Horizontal Tail

a0_h = 0.1071233*180/np.pi   #1/rad4.49
S_h = 18.25    #m^2
CL_h_design = -0.0821
M_crit_a_h = tl.bisection(tl.M_cr_calc_h,0.2,0.99,100)
# print(M_crit_a_h)
A_h =  3           # Design Choice
taper_h = 0.85      # Design Choice
sweep_le_h = np.radians(35)     # deg

# Vertical Tail

a0_v = 6.195564526     #1/rad
S_v = 16.9     # m^2
CL_v_design = 0.3986
M_crit_a_v = tl.bisection(tl.M_cr_calc_v,0.2,0.99,100)
# print(M_crit_a_v)
A_v =   1.1         # Design Choice
taper_v = 0.8     # Design Choice
sweep_le_v = np.radians(30)    # deg


"""Computations"""

# Horizontal Tail



b_h = np.sqrt(A_h*S_h)
MAC_h = S_h/b_h
Cr_h = (2*S_h)/(b_h*(1+taper_h))
# Cr_h = (3/2)*MAC_h*(1+taper_h)/(1+taper_h+(taper_h**2))
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
MAC_v = S_v/b_v
Cr_v = (2*S_v)/(b_v*(1+taper_v))
# Cr_v = (3/2)*MAC_v*(1+taper_v)/(1+taper_v+(taper_v**2))
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
print("MAC =", MAC_h)
print("M_crit =",M_crit_h, val_h)
print("Incidence angle =", np.rad2deg(CL_h_design/CL_h_alpha), "\n \n")


print("Vertical Tail \n \n")
print("Span =", b_v)
print("CL_alpha =", CL_v_alpha)
print("Cr =", Cr_v)
print("Ct =",Ct_v)
print("MAC =", MAC_v)
print("M_crit =",M_crit_v, val_v)
print("Required angle at OEI", np.rad2deg(CL_v_design/CL_v_alpha))

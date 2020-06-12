# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:19:31 2020

@author: keike
"""
import numpy as np

#Constants
K_c = 1.4

from Loading_diagram import xlemac, cg_aft, lmac, cg_frw, l_t

#Inputs
V_h = 1.1       # Tail volume assumed, M. Adraey [m^3]
D_f = 4.2       # Diameter fuselage [m]
S = 139.1       # Wing area [m^2]
AR = 17
lamda_LE = 40.63/360*2*np.pi
Cm_af = -0.1125

dH_mostfrw = (-(xlemac-0.3)+cg_frw)/lmac
dH_mostaft = (-(xlemac-0.3)+cg_aft)/lmac #most positive (X_cg-X_ac)/C 
dh = dH_mostaft
CL_cruise = 5.66E-01


S_h = 0.25*S
CMac = -0.261601
CL_h = (CMac + CL_cruise*dh)/(0.25*l_t/lmac)


downwash = 0
CL_a_wf = 4.1
CL_alpha_h = 4

CM_a = CL_a_wf*(dH_mostaft)- CL_alpha_h*S_h/S*((l_t-cg_aft)/lmac)*(1-downwash)
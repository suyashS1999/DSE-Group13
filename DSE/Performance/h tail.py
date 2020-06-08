# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:19:31 2020

@author: keike
"""
import numpy as np

#Constants
K_c = 1.4

from Loading_diagram import xlemac, cg_aft, lmac

#Inputs
V_h = 1.1       # Tail volume assumed, M. Adraey [m^3]
D_f = 4.2       # Diameter fuselage [m]
S = 139.1       # Wing area [m^2]
AR = 17
lamda_LE = 40.63/360*2*np.pi
Cm_af = -0.1125
dh = -0.1*lmac
CL_cruise = 5.66E-01

#Tail volume from M. Sadraey
#CL_H = (Cm_ac+CL_w*0.2*MAC)/(V_H)
l_h = K_c*(4*lmac*S*V_h/(np.pi*D_f))**(1/2)
S_h = V_h*lmac*S/l_h
CMac = Cm_af*(AR*(np.cos(lamda_LE))**2/(AR+2*np.cos(lamda_LE)))

CL_h = (CMac + CL_cruise*dh)/V_h


dH_mostaft = (-(xlemac+0.25*lmac)+cg_aft)/lmac #most positive (X_cg-X_ac)/C 
downwash = 0.3
l_t = 16
CL_a_wf = 4.1
CL_alpha_h = 4

CM_a = CL_a_wf*(dH_mostaft)- CL_alpha_h*S_h/S*((l_t-cg_aft)/lmac)*(1-downwash)
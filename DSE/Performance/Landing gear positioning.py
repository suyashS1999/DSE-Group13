# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 16:48:51 2020

@author: keike
"""

import numpy as np

#%% Inputs
CL_a = 4.0          #CL_alpha (more critical out of in- outboard)
CL_max_clean = 1.35 
CL_cruise = 0.566   
p = 0.15            #Correction factor Torenbeek (0.15 - 0.20 range)
CL_max_TO = 2.0     
dtheta_dt = 4       # Torenbeek correction factor (rotation right after LOF) [deg/sec]
V_min = 65.6        # V_min with CL_max_TO
V_LOF = 1.05*V_min  # Lift off velocity [m/s]
X = 2.1             # Length of the LG struts [m] (from bottom of fuselage to bottom of tires)
H_fuse = 4.2
Y_LG = 3.3
Z_cg = X + H_fuse*0.5

#Position of landing gears from nose
x_mlg = 21.8    # main [m]
x_nlg = 3.2     # nose [m]
from Loading_diagram import cg_frw, cg_aft #import most frw and aft cg positions

#%% Functions
def th_LOF(theta):
    """
    This function computes lift off pitch angle of given initial estimate
    Start at around theta = 15 deg
    """
    l = X/np.sin(theta*np.pi/180)
    a_LOF = 1/CL_a *(CL_max_clean - CL_cruise - p*CL_max_TO) #AoA during lift off [rad]
    theta_LOF = a_LOF*180/np.pi + dtheta_dt*(2*l/V_LOF+np.sqrt(l*CL_max_TO/9.81/CL_a))
    return theta_LOF-theta

def iterate(theta):
    """
    This function preforms iterations to find the lift off pitch angle
    Outputs the iterated theta [deg]
    """
    dth = th_LOF(theta)
    while abs(dth) > 0.01:
        theta = theta + dth
        dth = th_LOF(theta)
        print("iterating")
    return theta

def nose_load():
    frac_min = 1/((cg_frw - x_nlg)/(x_mlg - cg_frw)+1)
    frac_max = 1/((cg_aft - x_nlg)/(x_mlg - cg_aft)+1)
    if frac_min > 0.08 and frac_min < 0.15 and frac_max > 0.08 and frac_max < 0.15:
        print("Nose load within 8-15% range", round(frac_min,3), round(frac_max,3))
    else:
        print("nose load out of range, need redesign", round(frac_min,3), round(frac_max,3))
    return frac_min, frac_max

def lateral_LG():
    c = (cg_aft-x_nlg)*np.sin(Y_LG/(x_mlg-x_nlg))
    
    phi  = round(np.arctan(Z_cg/c)*180/np.pi, 1)
    if phi > 55:
        print("phi needs to be lower than 55 deg! phi is now", phi, "deg")
    else:
        print("phi is", phi, "deg")
    return phi
nose_load()
lateral_LG()
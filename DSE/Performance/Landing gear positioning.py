# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 16:48:51 2020

@author: keike
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Inputs
CL_a = 5.03          #CL_alpha (more critical out of in- outboard)
CL_max_clean = 1.736 
CL_cruise = 0.5101828488
MTOW =    79584.74877
MLW = 79584.74877*0.85
CL_max_TO = 2.1     
dtheta_dt = 4       # Torenbeek correction factor (rotation right after LOF) [deg/sec]
V_min = np.sqrt(2*MTOW*9.81/(1.225*150.7*CL_max_TO))
#V_min = 65.6        # V_min with CL_max_TO
V_LOF = 1.05*V_min  # Lift off velocity [m/s]
X = 2             # Length of the LG struts [m] (from bottom of fuselage to bottom of tires)
H_fuse = 4.37
Y_LG = 4.1
Z_cg = X + H_fuse*0.5
W_S = 4958
p = 0.15            #Correction factor Torenbeek (0.15 - 0.20 range)

#Standard positioning --> MLG 20.6, NLG 1.5, X = 3.1, Y = 4.2

#Position of landing gears from nose
x_mlg = 24.2    # main [m]
x_nlg = 2.3     # nose [m]
from Loading_diagram import cg_frw, cg_aft #import most frw and aft cg positions
#from Class_II_Torenbeek import MTOW
#cg_frw = cg_frw + 3
#cg_aft = cg_aft + 3
MTOW = 78953.35625
#%% Functions
def th_LOF(theta):
    """
    This function computes lift off pitch angle of given initial estimate
    Start at around theta = 15 deg
    """
    l = (X+3-0.45)/np.sin(theta*np.pi/180)
    a_LOF = 1/CL_a *(CL_max_clean - CL_cruise - p*CL_max_TO) #AoA during lift off [rad]
    theta_LOF = a_LOF*180/np.pi + dtheta_dt*(2*2.5/V_LOF+np.sqrt(l*CL_max_TO/9.81/CL_a))
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

def static_load():
    l_n = cg_aft-x_nlg
    l_m = x_mlg-cg_aft
    P_m = l_n/(l_m+l_n)*MTOW
    return P_m

def illust(theta):
    fig = plt.figure(figsize = (16, 4))
    plt.ylim(-5, 6)
    plt.plot([0, 4.2, 36.2, 42.6, 32, 0], [0, 4.37, 4.37, 3, 0, 0])
    plt.plot([42.6, 42.6+1.2, 42.6+1.2, 42.6, 42.6], [3+0.45, 3+0.45, 3-0.45, 3-0.45, 3+0.45])
    plt.scatter([cg_frw, cg_aft], [2.1, 2.1])
    plt.plot([cg_aft, 42], [2.1, 2.1+(42-cg_aft)/-np.tan(theta/180*np.pi)])
    plt.plot([42.6+1.2, 10], [3-0.45, 3-0.45+np.tan(theta/180*np.pi)*(10- 42.6+1.2) ])
    plt.scatter([x_mlg, x_nlg],[-X, -X])
    plt.scatter(10,-X)
    plt.grid()
    return

#def min_leg():
#    w = 0.9*(5262)**0.25
#    lamda = 2 #or 2.5
#    eta_s = 0.6 #or 0.65
#    eta_t = 0.47
#    
#    Dt = 460 #or 49
#    bt = 160 #or 19
#    S_t = lamda*MTOW*0.92/4/(p*np.sqrt(Dt*bt))
#    S = 1/eta_s*(w**2/(1.84*9.81*lamda)-eta_t*S_t)
#    return S

w = 0.9*(4790/9.81)**0.25
lamda = 2.0 #or 2.5
eta_s = 0.8 #or oleo-pneumatic absorbers
eta_t = 0.47
Dt = 46*0.0254 #or 49
bt = 17*0.0254 #or 19
pres = 1.379e+6
S_t = lamda*MLW*9.81*0.92/4/(pres*np.sqrt(Dt*bt))
S = 1/eta_s*(w**2/(1.84*9.81*lamda)-eta_t*S_t)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:15:46 2020

@author: youssef
"""

hscr = 50    # screen height [ft]
G    = 9.81
mu_free  = 0.03  # friction coefficient free rolling (0.03-0.05)
mu_brake =  0.3  # friction coefficient braking (0.3-0.5)

# estimate drag contribution due to flap deflection

Rf = 0.25   # flap chord/wing chord
Swet_to_Sref = 0.65  # ratio between wetted area and wing surface area
Sflap_to_Sref =Rf*Swet_to_Sref
defl = 50 #degrees

delta1 = 179.32*Rf**4 - 111.6*Rf**3 + 28.929*Rf**2 + 2.3705*Rf - 0.0089
delta2 = (-3.9877*10**-12)*defl**6 + (1.1685*10**-9)*defl**5 + (-1.2846*10**-7)*defl**4 + (6.1742*10**-6)*defl**3 +(-9.89444*10**-5)*defl**2 + (6.8324*10**-4)*defl**-4 + (-3.892*10**-4)

deltaCDflap = delta1*delta2*Sflap_to_Sref




# APPROACH DISTANCE


# FLARE DISTANCE


# FREE-ROLL DISTANCE

# BRAKING DISTANCE


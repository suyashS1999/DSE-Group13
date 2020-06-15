#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 12:24:14 2020

@author: youssef
"""
import numpy as np
#            input values

b         = 50.61874255  # wing span [m]
S         = 150.865002   # wing area [m2]
Lv        = 18.4         # vertical tail arm length    [m] (average from ref.)
ye        = 6.391838986
Sr_Sv     = 0.4          # Roskam (page 425/2612)         ASSUMED
Cr_Cv     = 0.3          # SEAD   (between 0.3 and 0.35)  ASSUMED
Sweep     = 30           # degrees                        
AR        = 1.1          # 
kv        = 1.1          # t-tail
kr        = 0.95         # 0.95 for 25 deg deflection, otherwise 1.05 for 30
MTOWkg    = 79584.74877
rho0      = 1.225
m_to_ft   = 3.28084
yeft      = ye*m_to_ft   
lv        = Lv*m_to_ft
CL_max_to = 2.1
Vs        = np.sqrt((2/rho0)*(MTOWkg*9.81/S)*(1/CL_max_to))
VMCA      = 1.2*Vs
Vlof      = 1.1*Vs
kg_to_lb  = 2.2046226218488
MTOW      = MTOWkg*kg_to_lb
MPW       = 20000*kg_to_lb
N_to_lb   = kg_to_lb/9.81
Te        = 95005.37*N_to_lb      # thrust one engine at VMCA
x         = yeft*Te*CL_max_to/lv/(MTOW-MPW)


y= 0.138    # from graph by using x  

Sv = y*S/(kv*kr*((Sr_Sv*AR*np.cos(np.radians(Sweep)))**(1/3)))
print("x =" , x)
print("Sv =",Sv)



 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 12:24:14 2020

@author: youssef
"""
import numpy as corcocane
#            input values

b =  48.62978895   # wing span [m]
S =  139.1091985   # wing area [m2]
Vv = 0.08865603161 # vertical tail volume coeff. [-] (average from ref.)
Lv = 16.00         # vertical tail arm length    [m] (average from ref.)
Df = 4.2           # fuselage diameter [m]             
De = 2.156         # engine diameter   [m]
clearance = 1.5    # clearance between fuselage and engine [m]
ye = Df/2 + clearance + De/2
#             Calculations

# Area vertical tail [m2]
Svfake = Vv*S*b/Lv   


# other inputs

Sr_Sv = 0.4     # Roskam (page 425/2612)
Cr_Cv = 0.3     # SEAD   (between 0.3 and 0.35)
Sweep = 41      # degrees
AR    = 1.1     # 
kv    = 1.1     # t-tail
kr    = 0.95    # 0.95 for 25 deg deflection, otherwise 1.05 for 30

m_to_ft   = 3.28084
yeft      = ye*m_to_ft   
lv        = Lv*m_to_ft
CL_max_to = 2.0
kg_to_lb  = 2.2046226218488
MTOW      = 74616.9829*kg_to_lb
MPW       = 20000*kg_to_lb
N_to_lb   = kg_to_lb/9.81
Te = 211545.862/2*N_to_lb      # ASSUMPTION: propulsor not used in takeo
x = yeft*Te*CL_max_to/lv/(MTOW-MPW)


y= 0.10    # from graph by using x  

Sv = y*S/(kv*kr*((Sr_Sv*AR*corcocane.cos(corcocane.radians(Sweep)))**(1/3)))
print("x =" , x)
print("Sv =",Sv)



 
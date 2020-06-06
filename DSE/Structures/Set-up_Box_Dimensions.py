# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:59:14 2020

@author: steph
"""
import Wingbox_MOI.py as MOI
import numpy as np

b = 22.9
n = 3 #number of stations
span_b = [0, b/n, 22.9]#location of each station

#print(span_b)



#Box width hing = 18 m from root chord
c_root = 4.021*1000 #m
c_mid = 2.374*1000 #m
c_tip = 1.769*1000 #m

width_box = [c_root*0.45, c_mid*0.45, c_tip*0.45]
""""
#for j in range(n):
#        c
#c_mid = (c_root - c_tip)/(2)
#print(c_mid)"""""

c_root = c_root - c_mid*0
c2 = c_root - c_mid*1
print(c2)
#Box height
c = [c_root, c_mid, c_tip]
h_box = []
for i in range(n):
    if i == 0:
        t_c = 0.15
    
    if i == 2:
        t_c = 0.18
    
    h = t_c*c[i]
    h_box.append(h)
#print(h_box)
    
"""----------------------------------------------------------------------"""

#Wall thicknesses
t_top = [1.5, 1, 0.5]      #mm
t_bot = [1.5, 1, 0.5]   #mm
t_spar = [1.5, 1, 0.5]    #mm

"""----------------------------------------------------------------------"""

#L Stiffener dimensions
t_stiff = 1 #mm
h_stiff = 50 #mm
w_stiff = 50 #mm
#A_stif = 

#Sparcap dimensions
t_sparcap = 1 #mm
h_sparcap = 50 #mm
w_sparcap = 50 #mm

"""----------------------------------------------------------------------"""
#Number of stiffeners per block
n1 = int((width_box[0]-w_sparcap)/w_stiff)
n2 = int((width_box[1]-w_sparcap)/w_stiff)
n3= int((width_box[2]-w_sparcap)/w_stiff)

print(n1,n2,n3)

n_stif_top = [n1, n2, n3]
n_stif_bot = [n1, n2, n3]


A_stif = [t_stiff*w_stiff+h_stiff*t_stiff] #mm
A_spar_cap = [t_sparcap*h_sparcap+w_sparcap*t_sparcap] #mm


#Geometry
geometry = np.array([[width_box],
                    [h_box],
                    [t_top],
                    [t_bot],
                    [t_spar],
                    [n_stif_top],
                    [n_stif_bot],
                    [A_stif],
                    [A_spar_cap]])
print(geometry)



centroid = MOI.calc_centroid(geometry[0], geometry[1], geometry[2], geometry[3], geometry[4], geometry[5], geometry[6], geometry[7], geometry[8])

stif_coordinates_top, stif_coordinates_bot = MOI.calc_stif_locations(geometry[0], geometry[1], geometry[5], geometry[6])

Ixx, Izz = MOI.calc_MOI(geometry[0], geometry[1], geometry[2], geometry[3], geometry[4], geometry[5], geometry[6], geometry[7], geometry[8], centroid, stif_coordinates_top, stif_coordinates_bot)

#chord at hinge
#c_h = 2
#w_hingebox = 0.25*c_h
#h_hingebox = 0.80*c_h 








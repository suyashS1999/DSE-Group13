# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 13:48:53 2020

@author: keike
"""

import numpy as np
import matplotlib.pyplot as plt

#General position information
lmac = 1.501
xlemac = 4 + (32.05*0.5) - 0.25*lmac          #quarter MAC at 50% of cabin (SEAD rough approximation)
x_cargo_frw = 7.07
x_cargo_aft = 37.7
x_cg_fuel = xlemac

#CG locations
l_fuse = 42.6
x_cg_wing = xlemac+lmac*0.25 #?
x_cg_tail = 5.75
x_cg_body = 20.08
x_cg_undercarriage_main = 0.6*l_fuse
x_cg_undercarriage_nose = 0.1*l_fuse
x_cg_surfcont = x_cg_wing
x_cg_nacelle = xlemac
x_cg_engine = x_cg_nacelle
x_cg_equip = x_cg_body

#Determine CG of the OEW w.r.t. nose
from Class_II_Torenbeek import *
x_cg_oew = (W_wing*x_cg_wing + W_tail*x_cg_tail + W_body*x_cg_body + W_undercarriage_main*x_cg_undercarriage_main + W_undercarriage_nose*x_cg_undercarriage_nose + W_surfcont*x_cg_surfcont + W_nacelle*x_cg_nacelle + W_engine*x_cg_engine + W_equip*x_cg_equip)/OEW

#Passenger arrangement (input value)
total_pax = 192
window_row = 2
middle_row = 2
aisle_row = 2
seat_pitch = 812.8e-3 #[m] 
x_seat_mostfrw = 7 #[m] Assume all seat at the middle
pax_weight = 95*194/192 #[kg]

#General weight information (input value)
#OEW = 42145.77058 #[kg] #Including pilots and crew etc!
PL = 20000
m_cargo_frw = (PL - total_pax*pax_weight)*4/7
m_cargo_aft = (PL - total_pax*pax_weight)*3/7
m_fuel = 12471.21232


#Determine passenger number per row
pax_per_row = total_pax/(window_row + middle_row + aisle_row)

if pax_per_row == int(pax_per_row):
    pax_per_row = int(pax_per_row)
else:
    print("number of passenger per row is not equal!")
    exit()

def cargo_component_cg(m_frw, m_aft, x_frw, x_aft):
    frw_aft = [[],[]]
    aft_frw = [[],[]]
    frw_aft[0] = [m_frw, m_aft]
    frw_aft[1] = [x_frw, x_aft]
    aft_frw[0] = [m_aft, m_frw]
    aft_frw[1] = [x_aft, x_frw]
    return frw_aft, aft_frw

def passenger_component_cg(number_of_rows):
    frw_aft = [[],[]]
    aft_frw = [[],[]]
    m = np.ones(pax_per_row)*pax_weight*number_of_rows
    for i in range(pax_per_row):
        frw_aft[1].append(x_seat_mostfrw + i*seat_pitch)
    aft_frw[1] = sorted(frw_aft[1], reverse=True)
    frw_aft[0] = m
    aft_frw[0] = m
    return frw_aft, aft_frw

#Potato diagram function definition
"""
Input:  frw_aft, aft_frw (each consist of two lists, one list of mass of each component, and one list of x_loc of each component)
        m_initial, x_cg_initial (initial mass and x_cg for building up the potato diagram)
Output: x_cg_frw, x_cg_aft (the shitf of x_cg by adding the component one by one. frw: forward to backwards, aft: other way around)
        mlst (accumulated mass, one list assuming all component has equal mass)
             (which is not the case for the cargo! so make the second mlst by yourself which should be easy enough)
"""
def potato(frw_aft, aft_frw, m_initial, x_cg_initial): 
    m_current = m_initial
    mom_current = m_initial*x_cg_initial
    x_cg_frw = [x_cg_initial]
    for i in range(len(frw_aft[0])):        #Find xcg frw-->aft
        m_current = m_current + frw_aft[0][i]
        mom_current = mom_current + frw_aft[0][i]*frw_aft[1][i]
        x_cg_current = mom_current/m_current
        x_cg_frw.append(x_cg_current)
    m_current = m_initial
    mom_current = m_initial*x_cg_initial
    x_cg_aft = [x_cg_initial]
    for i in range(len(aft_frw[0])):        #Find xcg aft-->frw
        m_current = m_current + aft_frw[0][i]
        mom_current = mom_current + aft_frw[0][i]*aft_frw[1][i]
        x_cg_current = mom_current/m_current
        x_cg_aft.append(x_cg_current)    
    mlst = [m_initial]
    for i in range(len(frw_aft[0])):
        mlst.append(mlst[i]+frw_aft[0][i])
    return x_cg_frw, x_cg_aft, mlst #For cargo, mlst frw and aft are different it is just reverse of each other


###Main part --- plotting everything###
cargo = cargo_component_cg(m_cargo_frw, m_cargo_aft, x_cargo_frw, x_cargo_aft)
pass_window = passenger_component_cg(window_row)
pass_aisle  = passenger_component_cg(aisle_row)
pass_middle = passenger_component_cg(middle_row)

#Plotting the position with x from nose
cargo_potato  = potato(cargo[0], cargo[1], OEW, x_cg_oew)
window_potato = potato(pass_window[0], pass_window[1], cargo_potato[2][-1],  cargo_potato[0][-1])
aisle_potato  = potato(pass_aisle[0],  pass_aisle[1],  window_potato[2][-1], window_potato[0][-1])
middle_potato = potato(pass_middle[0], pass_middle[1], aisle_potato[2][-1],  aisle_potato[0][-1])
#Compute the mass change for cargo, it is asymmetrical
cargo_potato_mass = [cargo_potato[2], [cargo_potato[2][0], cargo_potato[2][0]+m_cargo_frw, cargo_potato[2][2]]]
#Plotting fuel manually
fuel_mass_shift = [middle_potato[2][-1], middle_potato[2][-1]+m_fuel]
x_cg_shift_fuel = [middle_potato[0][-1], (middle_potato[0][-1]*middle_potato[2][-1]+x_cg_fuel*m_fuel)/(middle_potato[2][-1]+m_fuel)]


#Convert into percentage of the mac from lemac
def pcmac(lst):
    lst = (np.array(lst)-xlemac)/lmac
    return lst

cargo_frw = pcmac(cargo_potato[0])
cargo_aft = pcmac(cargo_potato[1])
window_frw = pcmac(window_potato[0])
window_aft = pcmac(window_potato[1])
aisle_frw = pcmac(aisle_potato[0])
aisle_aft = pcmac(aisle_potato[1])
middle_frw = pcmac(middle_potato[0])
middle_aft = pcmac(middle_potato[1])
x_cg_shift_fuel = pcmac(x_cg_shift_fuel)

plt.plot(cargo_frw, cargo_potato_mass[0], label = "Cargo", color = "purple")
plt.plot(cargo_aft, cargo_potato_mass[1],  color = "purple")
plt.plot(window_frw, window_potato[2], label = "Window", color = "green")
plt.plot(window_aft, window_potato[2], color = "green")
plt.plot(aisle_frw, aisle_potato[2], label = "Aisle", color = "blue")
plt.plot(aisle_aft, aisle_potato[2], color = "blue")
plt.plot(middle_frw, middle_potato[2], label = "Middle", color = "red")
plt.plot(middle_aft, middle_potato[2], color = "red")
plt.plot(x_cg_shift_fuel, fuel_mass_shift, label = "Fuel", color = "orange")
plt.title("Loading diagram of Aircraft")
plt.legend()
plt.xlabel("x_cg [mac]")
plt.ylabel("Mass [kg]")
plt.show()
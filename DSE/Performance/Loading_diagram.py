# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 13:48:53 2020

@author: keike
"""

import numpy as np
import matplotlib.pyplot as plt

#%% ---------------------- Inputs ----------------------


#Aircraft dimensions
lmac = 1.544                          # MAC [m]
l_fuse = 42.6                           # Fuselage length [m]
xlemac = 0.46*l_fuse    # X position of LEMAC from nose (LEMAC at 50% of cabin)
x_cg_wing = xlemac + 0.25*lmac            # wing cg from nose, assumed same as average aerodynamic centre
#l_t = 42.6-6*0.75-1-x_cg_wing+1              # Tail moment arm [m], distance between 1/4 chords tail to xac (for now V_tail)
W_wing = 9235.824024                        # Wing and strut mass combined
W_engine = 983.4235289 + 4341.959989     # Engine + generator + nacelle
x_cg_engine = xlemac - 6.167                # assumed to be same as pylons

#CG locations in X direction from the nose (Assumptions/Justification in brackets)
x_cargo_frw = 7                        # forward cargo [m] (from concept 3 midterm side view) 
x_cargo_aft = 27                      # aft cargo [m] (from concept 3 midterm side view)
x_cg_fuel = xlemac                      # fuel [m] (approximated to be at the LEMAC, forward of approximate CG of wing since there will be no fuel at the tips)

#From Class II weight estimation google sheet
x_cg_oew_lesswing = 21.22121347 #[m] from nose
OEW_lesswing = 26344.20657
OEW = OEW_lesswing+W_wing+W_engine #[kg]
x_cg_oew = (x_cg_oew_lesswing*OEW_lesswing + x_cg_wing*W_wing + x_cg_engine*W_engine)/OEW
#Verify OEW matches

#Commented out, because this is calculated now in excel
#x_cg_wing = xlemac-10.58+9            # wing, as function of MAC [m] (For now assumed to be quarter chord MAC)
#x_cg_tail = l_t+xlemac               # tail [m] (cg wing + tail arm)
#x_cg_body = l_fuse*0.5                  # body/fuselage [m] (midpoint of fuselage length)
#x_cg_undercarriage_main = 0.53*l_fuse   # main landing gear [m] (relative to fuselage length of A320)
#x_cg_undercarriage_nose = 0.095*l_fuse  # nose landing gear [m] (relative to fuselage length of A320)
#x_cg_surfcont = x_cg_wing               # surface control [m] (assumed to be same as wing cg)
#x_cg_nacelle = xlemac-10.58+9                   # nacelle + pylons [m] (assumed to be at LEMAC)
#x_cg_engine = xlemac-10.58+9                    # engines [m] (assumed to be at LEMAC)
#x_cg_equip = x_cg_body                  # equipment e.g. APU and all that stuff [m] (assumed to be equally spread along the body)
#x_cg_furnish = x_cg_body
#x_cg_BLI = l_fuse

#Passenger arrangement (input value)
total_pax = 192                         # Number of passengers [-] (ITS 194 BUT MODIFIED SO THAT ROWS FIT NICELY)
window_row = 2                          # number of rows [-] (all from mid-term)
middle_row = 2
aisle_row = 2
total_row = window_row+middle_row+aisle_row
pax_per_row = int(total_pax/total_row)
seat_pitch = 812.8e-3*1.1                  # seating pitch [m] (from mid-term)
l_cabin = 30.25                         # length of cabin [m] (from mid-term)
l_cockpit = 4                           # length of cockpit [m] (from mid-term)
x_seat_mostfrw = l_cockpit+(l_cabin-(seat_pitch*total_pax/(total_row)))/2 # Position of most forward seat from nose [m] (Assume all seat at the middle of cabin with the seating pitch)
pax_weight = 75*194/192                 # weight per pax [kg] (assumed 95, but corrected to 192 pax from 194 pax)

# Aircraft weights
PL = 20000                                      # Payload mass [kg] (from req.)
m_cargo_frw = (PL - total_pax*pax_weight)*4/7   # Frw cargo [kg]
m_cargo_aft = (PL - total_pax*pax_weight)*3/7   # Aft cargo [kg]
m_fuel = 13765.99377                            # Fuel weight [kg]
#Verified that PL mass indeed equals PAX+cargofrw+cargoaft
#%% ---------------------- Main ----------------------
#Compile CGs
#from Class_II_Torenbeek import *
#x_cg_oew = (W_wing*x_cg_wing + W_tail*x_cg_tail + W_body*x_cg_body + W_undercarriage_main*x_cg_undercarriage_main + W_undercarriage_nose*x_cg_undercarriage_nose + W_surfcont*x_cg_surfcont + W_nacelle*x_cg_nacelle + W_engine*x_cg_engine + W_equip*x_cg_equip + W_furnish*x_cg_furnish + W_BLI*x_cg_BLI)/OEW

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

#Convert into percentage of the mac from lemac
def pcmac(lst):
    lst = (np.array(lst)-xlemac)/lmac
    return lst


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

cargo_frw = pcmac(cargo_potato[0])
cargo_aft = pcmac(cargo_potato[1])
window_frw = pcmac(window_potato[0])
window_aft = pcmac(window_potato[1])
aisle_frw = pcmac(aisle_potato[0])
aisle_aft = pcmac(aisle_potato[1])
middle_frw = pcmac(middle_potato[0])
middle_aft = pcmac(middle_potato[1])
x_cg_shift_fuel = pcmac(x_cg_shift_fuel)

def loading_diagram():
    plt.plot(cargo_frw, cargo_potato_mass[0], label = "Cargo", color = "purple")
    plt.plot(cargo_aft, cargo_potato_mass[1],  color = "purple")
    plt.plot(window_frw, window_potato[2], label = "Window", color = "green")
    plt.plot(window_aft, window_potato[2], color = "green")
    plt.plot(aisle_frw, aisle_potato[2], label = "Aisle", color = "blue")
    plt.plot(aisle_aft, aisle_potato[2], color = "blue")
    plt.plot(middle_frw, middle_potato[2], label = "Middle", color = "red")
    plt.plot(middle_aft, middle_potato[2], color = "red")
    plt.plot(x_cg_shift_fuel, fuel_mass_shift, label = "Fuel", color = "orange")
#    plt.title("Loading diagram of Aircraft")
    plt.legend(fontsize = 14, loc = 'upper left')
    plt.xlabel("CG location of X direction [MAC]", fontsize = 14)
    plt.ylabel("Mass [kg]", fontsize = 14)
#    plt.xlim(-0.55, 0.45)
    plt.ylim(OEW, 82000)
    plt.grid()
    plt.show()
    return 
    

if __name__ == "__main__":
    loading_diagram()

cg_frw_mac = np.min([np.min(cargo_frw), np.min(window_frw), np.min(aisle_frw), np.min(middle_frw), np.min(x_cg_shift_fuel)])
cg_aft_mac = np.max([np.max(cargo_aft), np.max(window_aft), np.max(aisle_aft), np.max(middle_aft), np.max(x_cg_shift_fuel)])
cg_frw = xlemac + cg_frw_mac*lmac
cg_aft = xlemac + cg_aft_mac*lmac
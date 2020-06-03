#%% ---------------------- About this code ----------------------
# The coordinate system is such that the lengths from nose gets more positive as you go more aft


#%% ---------------------- Imports ----------------------


#%% ---------------------- Constants ----------------------


#%% ---------------------- Inputs ----------------------
MAC = 5         # (DUMMY) Mean aerodynamic chord [m]
LE_MAC = 15     # (DUMMY) Position of leading edge of MAC from nose [m]
CG_wing = 3     # (DUMMY) Position of CG of the main wing w.r.t. LEMAC [m]
CG_engine = 2   # (DUMMY) Position of CG of the engines w.r.t. LEMAC [m]

#The weights [kg] and location of CG from the nose [m]
wing  = [500, LE_MAC+CG_wing]           # (DUMMY) main wing
strut = [200, 17]
nacelle_pylon = [200, LE_MAC+CG_engine] # (DUMMY) nacelle and pylons
v_tail = [300, 40]
h_tail = [250, 41]
fuselage = [5000, 20]
#%% ---------------------- About this code ----------------------
# This program generates the flight envelope of the aircraft
# Inputs are:
# Outputs are:


#%% ---------------------- Imports ----------------------
import numpy as np
import matplotlib.pyplot as plt

#%% ---------------------- Constants ----------------------
g = 9.81                        # Gravitational acceleration [m/s^2]
M_cruise = 0.78                 # Mach number in cruise [-]
h_cruise = 11000                # Altitude in cruise [m]
rho0 = 1.225                    # Sea level air Density [kg/m^3]

#%% ---------------------- Inputs ----------------------
MTOW = 74616.9829*g             # MTOW [N]
S_ref =  142.520611             # Reference wing area [m^2]
CL_max_clean = 1                # DUMMY Max CL with clean configuration
n1 = 2.5                        # Positive limit maneuvering load factor [-]
#%% ---------------------- Functions ----------------------
#Copy pasted for now, cannot import .py files inside other folders
def ISA_trop(h):
	""" This function computes the atmospheric properties 
		within the troposphere
	Input:
		h = altitude [m]
	Output:
		T = Temperature [K]
		p = Pressure [Pa]
		rho = Densituy [kg/m^3]
		a = spped of sound [m/s]
	"""
	T = 288.15 - 0.0065*h;
	p = 101325*(T/288.15)**(-g/(-0.0065*287));
	rho = 1.225*(T/288.15)**(-g/(-0.0065*287) - 1);
	a = sqrt(1.4*287*T);
	return T, p, rho, a;


#%% ---------------------- Main ----------------------
T, p, rho, a = ISA_trop(h)
W_cruise = MTOW;                                            # For now I assume W is MTOW during cruise
V_stall = np.sqrt(W_cruise/S_ref*2/rho/CL_max_clean)        # Stall speed [m/s]
V_maneuver = V_stall*np.sqrt(n1)                            # Maneuver speed [m/s]
V_cruise = a*M_cruise                                       # Cruise speed as function of M [m/s]
V_dive = V_cruise/0.8                                      # Dive speed

#%% ---------------------- About this code ----------------------
# This program generates the flight envelope of the aircraft
# Inputs are:
# Outputs are:


#%% ---------------------- Imports ----------------------
import numpy as np
import matplotlib.pyplot as plt
from ConceptEvaluation.PayloadRangeDiagram import ISA_trop

#%% ---------------------- Constants ----------------------
g = 9.81                        # Gravitational acceleration [m/s^2]
M_cruise = 0.78                 # Mach number in cruise [-]
h_cruise = 11000                # Altitude in cruise [m]

#%% ---------------------- Inputs ----------------------
MTOW = 74616.9829                   # MTOW [kg]
S_ref =  142.520611                 # Reference wing area [m^2]
CL_max_clean = 1.5                  # DUMMY Max CL with clean configuration [-]
CL_min_clean = -0.8*CL_max_clean    # DUMMY negative Max CL with clean configuration [-]
CL_max_flap = 2.0                   # MAX CL with LANDING CONFIGURATION (need to be determined) [-]
n_max = 2.5                         # Positive limit maneuvering load factor [-]
n_min = -1                          # Negative limit maneuvering load factor [-]
n_flap = 2                          # Assumed (from 2nd year project guide)
#%% ---------------------- Functions ----------------------



#%% ---------------------- Main ----------------------
T, p, rho, a = ISA_trop(h_cruise)
T0, p0, rho0, a0 = ISA_trop(0)

W_cruise = MTOW*g;                                          # For now I assume W is MTOW during cruise
V_stall_clean = np.sqrt(W_cruise/S_ref*2/rho0/CL_max_clean) # Stall speed EAS [m/s]
V_stall_flap = np.sqrt(W_cruise/S_ref*2/rho0/CL_max_flap)   # DUMMY Stall speed EAS with flaps extended
V_maneuver = V_stall_clean*np.sqrt(n_max)                   # Maneuver speed EAS [m/s]
V_cruise = a0*M_cruise*np.sqrt(p/p0)                        # Cruise speed EAS [m/s]
V_dive = V_cruise/0.8                                       # Dive speed EAS [m/s]
V_neg_stall_clean = np.sqrt(W_cruise/S_ref*2/rho0/-CL_min_clean)
V_max_flap = 0.85*V_maneuver                                 # Max speed with flaps extended

#Maneuver envelope
fig = plt.figure(figsize = (12, 8))
plt.xlabel("Equivalent airspeed, $V_{EAS}$")
plt.ylabel("Load factor, n")
plt.axhline(y = 0, color = "black")
plt.xlim(0, V_dive + 10)

#Limiting velocities EAS
plt.vlines(x = V_dive, ymin = 0, ymax = n_max)
plt.text(V_dive, 0, '$V_D$', ha='left', va='bottom')
#V_cruise and V_dive linear 
plt.plot([V_cruise, V_dive], [n_min, 0], color = "black")
#Origin to n_max clean
V_positive_clean = np.linspace(0, V_maneuver, 100)
n_positive_clean = (V_positive_clean/V_stall_clean)**2
plt.plot(V_positive_clean, n_positive_clean, color = "black")
#Origin to n_min
n_negative_clean = np.linspace(0, n_min, 100)
V_negative_clean = np.sqrt(-n_negative_clean)*V_neg_stall_clean
plt.plot(V_negative_clean, n_negative_clean, color = "black")
#Limiting loads
plt.hlines(y = n_max, xmin = V_maneuver, xmax = V_dive)
plt.hlines(y = n_min, xmin = V_neg_stall_clean, xmax = V_cruise)

#Flaps extended
#Origin to n_max flaps
n_positive_flap = np.linspace(0, n_flap, 100)
V_positive_flap = np.sqrt(n_positive_flap/2)*V_stall_flap
plt.plot(V_positive_flap, n_positive_flap, "--", color = "black")
#Horizontal line
plt.hlines(y = n_flap, xmin = V_positive_flap[-1], xmax = V_max_flap)
#Vertical line
plt.vlines(x = V_max_flap, ymin = 0, ymax = n_flap, linestyles = "--")
plt.text(V_max_flap, 0, '$V_F$', ha='left', va='bottom')

#Dotted lines to help make it easier to see
plt.vlines(x = V_cruise, ymin = n_min, ymax = n_max, linestyles = "--")
plt.text(V_cruise, 0, '$V_C$', ha='left', va='bottom')
plt.hlines(y = 1, xmin = 0, xmax = V_dive, linestyles = "--")
plt.vlines(x = V_stall_clean, ymin = 0, ymax = 1, linestyles = "--")
plt.text(V_stall_clean, 0, '$V_{S1}$', ha='left', va='bottom')
plt.vlines(x = V_stall_flap, ymin = 0, ymax = 2, linestyles = "--")
plt.text(V_stall_flap, 0, '$V_{S0}$', ha='left', va='bottom')
plt.vlines(x = V_maneuver, ymin = 0, ymax = (V_maneuver/V_stall_clean)**2, linestyles = "--")
plt.text(V_maneuver, 0, '$V_A$', ha='left', va='bottom')
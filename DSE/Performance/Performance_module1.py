#%% ---------------------- About this code ----------------------
# This program generates the flight envelope of the aircraft
# It makes both maneuver and gust V-n diagrams in one graph
# Inputs are:
# Outputs are:


#%% ---------------------- Imports ----------------------
import numpy as np
import matplotlib.pyplot as plt
#from ConceptEvaluation.PayloadRangeDiagram import ISA_trop

#%% ---------------------- Constants ----------------------
g = 9.81                        # Gravitational acceleration [m/s^2]
M_cruise = 0.78                 # Mach number in cruise [-]
h_cruise = 11000                # Altitude in cruise [m]
h = h_cruise

#%% ---------------------- Inputs ----------------------
MTOW = 71311.11553                   # MTOW - fuel used in ascend [kg]
S_ref =  139.1                 # Reference wing area [m^2]
CL_max_clean = 1.35                # Max CL with clean configuration [-]
CL_min_clean = -CL_max_clean        # Assumed same magnitude, negative Max CL with clean configuration [-] (not so critical)
CL_max_flap = 2.4                   # MAX CL with LANDING CONFIGURATION [-]
n_max = 2.5                         # Positive limit maneuvering load factor [-]
n_min = -1                          # Negative limit maneuvering load factor [-]
n_flap = 2                          # Assumed (from 2nd year project guide)

CL_alpha = 5.267/np.sqrt(1-M_cruise**2)                    # CL alpha gradient of inboard (higher gradient more critical) [-/rad]
W_S = 5262                          # Wing loading [N/m^2]
c = 1.501                           # Mean aerodynamic chord [m]
#%% ---------------------- Functions ----------------------
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
	a = np.sqrt(1.4*287*T);
	return T, p, rho, a;

def relevant_V(W_cruise):
    """
    This function outputs all the relevant characteristic velocities
    These are stall speeds, maneuver speeds, cruise speeds, dive speeds
    All of these velocities are EAS!!!
    Input is he weight of the aircraft [N]
    """
    V_stall_clean = np.sqrt(W_cruise/S_ref*2/rho0/CL_max_clean)     # Stall speed EAS [m/s]
    V_stall_flap = np.sqrt(W_cruise/S_ref*2/rho0/CL_max_flap)       # Stall speed EAS with flaps extended
    V_maneuver = V_stall_clean*np.sqrt(n_max)                       # Maneuver speed EAS [m/s]
    V_neg_stall_clean = np.sqrt(W_cruise/S_ref*2/rho0/-CL_min_clean)# Limiting velocity when n = -1 [m/s]
                               # Cruise speed EAS [m/s]
                                              # Dive speed EAS [m/s]
    return (V_stall_flap, V_stall_clean, V_neg_stall_clean, V_maneuver)

def gust_n(V_EAS, U_EAS):
    """
    This function gives the limiting load factors due to gust with given velocity
    """
    delta_n = K_g*CL_alpha*rho0*U_EAS*V_EAS/(2*W_S)
    return (delta_n)
    
#%% ---------------------- Main ----------------------
T, p, rho, a = ISA_trop(h)
T0, p0, rho0, a0 = ISA_trop(0)
V_stall_flap, V_stall_clean, V_neg_stall_clean, V_maneuver = relevant_V(MTOW*g) # For now I assume W is MTOW during cruise
V_max_flap = 1.8*V_stall_flap                                 # Max speed with flaps extended EAS [m/s]
V_cruise = a0*M_cruise*np.sqrt(ISA_trop(h_cruise)[1]/p0) 
V_dive = V_cruise/0.8
#Maneuver envelope
fig = plt.figure(figsize = (12, 8))
plt.xlabel("Equivalent airspeed, $V_{EAS}$")
plt.ylabel("Load factor, n")
plt.axhline(y = 0, color = "black")
#plt.xlim(0, V_dive + 10)

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


###Gust envelope###
w = W_S*0.020885439                                         # W/S [lb/ft^2]
mu = 2*w/(rho*0.0019403203*c/0.3048*CL_alpha*g/0.3048)
K_g = 0.88*mu/(5.3+mu)                                      # Gust alleviation factor [-]
dU_refdH = (7.92-13.41)/(15240-4572)                        # U_ref relation with altititude higher than 4572 m [m/s /m]
U_ref = ((h_cruise-4572)*dU_refdH+13.41)/0.3048             # Gust reference velocity [ft/s] (CS25. 341(5))
V_alpha_max = V_stall_clean*(1+K_g*U_ref*V_cruise*1.94384*CL_alpha/(498*w))**(1/2) #[m/s]
U_alpha_max = K_g*51*0.3048                                 # Statistical gust velocity high alpha TAS [m/s]
U_cruise = K_g*36.5*0.3048                                  # Statistical gust velocity cruise TAS [m/s]
U_dive = K_g*18.5*0.3048                                    # Statistical gust velocity dive TAS [m/s]
U_gusts = np.array([U_alpha_max, U_cruise, U_dive])
V_gusts = np.array([V_alpha_max, V_cruise, V_dive])
n_gusts = gust_n(V_gusts, U_gusts)




plt.plot([0,V_alpha_max, V_cruise, V_dive], [1, 1+n_gusts[0], 1+n_gusts[1], 1+n_gusts[2]], color = "black")
plt.plot([0,V_alpha_max, V_cruise, V_dive], [1, 1-n_gusts[0], 1-n_gusts[1], 1-n_gusts[2]], color = "black")
plt.plot([0, V_alpha_max], [1, 1+n_gusts[0]], "--", color = "black")
plt.plot([0, V_cruise], [1, 1+n_gusts[1]], "--", color = "black")
plt.plot([0, V_dive], [1, 1+n_gusts[2]], "--", color = "black")
plt.plot([0, V_alpha_max], [1, 1-n_gusts[0]], "--", color = "black")
plt.plot([0, V_cruise], [1, 1-n_gusts[1]], "--", color = "black")
plt.plot([0, V_dive], [1, 1-n_gusts[2]], "--", color = "black")
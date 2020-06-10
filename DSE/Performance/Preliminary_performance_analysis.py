from ConceptEvaluation.PayloadRangeDiagram import *
from ConceptEvaluation.WingandThrustLoadingDiagrams import *

#%% ---------------------- Constants ----------------------
g = 9.81;										# Gravitaional acceleration [m/s^2]
rho0 = 1.225;									# Sea level air Density [kg/m^3]
M = 0.78;										# Cruise mach number [-]
#%% ------------------------ INPUTS ------------------------
MTOW = 76446.83536*g;							# Max Take off Weight [N]
MLW = 0.86*MTOW									# Max Landing Weight [N]
OEW = 42939.92655*g;							    # Operational empty weight [N]
Payload_max = 20*1000*g;						# Payload weight [N]
Max_Fuel_cap = 17.9*1000*g;						# Max fuel weight [N]
cj = 0.0000143;							        # Specific fuel consumption [kg/Ns]
L_D = 20.41115173;								# Lift to drag ratio [-]
h = 11000;										# Cruise altitude [m]
e = 0.7023214054;										# Ozwald efficiency [-]
landdis = 1800;									# Landing distance [m]
Cd0 = 0.0192;									# Zero lift drag coefficient [-]
k = 190;										# Take off parameter from Raymer and a/c database [retard units]
v_stall_landing = sqrt(landdis/0.5847);			# Stall speed calcuated emperically [m/s]
rho_airport = 1.225;							# air density at runway altitude [kg/m^3]
sigma = rho_airport/rho0;						# Density ratio [-]
N_engines = 2;									# Number of engines [-]
#c_v = 0.01199;	(now specified inside function)	# Climb gradient divied by velocity [-]
Cl_max = array([2.4, 2.5]);							# Cl max values for assessment [-]
Cl_max_takeoff = array([2.0, 2.1]);				 	# Cl_max for take off [-]
A = array([17]);								# Aspect ratio [-]
W_S_max = 7000;									# Max Wing Loading value, change this value if you want to change the range of wing loading values you want to assess [N/m^2]
n = 2.5;										# Maximum load factor [-]
pre_cruise_fuel_frac = 0.2507559723;					# Pre cruise fuel fraction [-]
post_cruise_fuel_frac = 0.08502789007;					# Post cruise fuel fraction [-]
reserve_fuel_frac = 0.1128409767;						# Reserve fuel fraction [-]
#%% ------------------------ Main ------------------------
T, p, rho_cruise, a = ISA_trop(h);
sigma_cruise = rho_cruise/rho0;			# Density ratio [-]
v_cruise = M*a;							# Cruise speed
PayloadRangeDiagram_JET(MTOW, OEW, Payload_max, reserve_fuel_frac, pre_cruise_fuel_frac, post_cruise_fuel_frac, Max_Fuel_cap, (g, M, a, cj, L_D));

fig = plt.figure(figsize = (10, 8));
_ = W_S_stall(v_stall_landing, Cl_max, MLW/MTOW, fig);
_, _ = W_S_takeoff(Cl_max_takeoff, k, sigma, W_S_max, fig);
_, _ = W_S_cruise(A, Cd0, rho_cruise, sigma_cruise, v_cruise, W_S_max, fig);
_ = W_S_climb_grad(Cd0, A, e, W_S_max, N_engines, fig);
#_, _ = W_S_maneuvering(n, Cd0, rho0, v_stall_landing, A, e, W_S_max, fig);

plt.grid(True);
plt.title("Thrust Loading vs Wing Loading (Concept 3)")
plt.axis([0, W_S_max, 0, 0.5]);
plt.legend(bbox_to_anchor=(0, 1), loc='upper left');
plt.xlabel("Wing Loading [N/m^2]");
plt.ylabel("Thrust Loading [-]");
plt.show();


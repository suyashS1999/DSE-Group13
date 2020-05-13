from PayloadRangeDiagram import*
from WingandThrustLoadingDiagrams import*

#%% ---------------------- Constants ----------------------
g = 9.81;										# Gravitaional acceleration [m/s^2]
rho0 = 1.225;									# Sea level air Density [kg/m^3]
M = 0.78;										# Cruise mach number [-]
#%% ------------------------ INPUTS ------------------------
MTOW = 83810.83645*g;							# Max Take off Weight [N]
MLW = 0.86*MTOW									# Max Landing Weight [N]
OEW = 46135.90302*g;							# Operational empty weight [N]
Payload_max = 20*1000*g;						# Payload weight [N]
Max_Fuel_cap = 21*1000*g;						# Max fuel weight [N]
cj = 0.00001416296856;							# Specific fuel consumption [kg/Ns]
L_D = 15.26974393;								# Lift to drag ratio [-]
h = 11000;										# Cruise altitude [m]
e = 0.8;										# Ozwald efficiency [-]
landdis = 1800;									# Landing distance [m]
Cd0 = 0.0192;									# Zero lift drag coefficient [-]
k = 190;										# Take off parameter from Raymer and a/c database [retard units]
v_stall_landing = sqrt(landdis/0.5847);			# Stall speed calcuated emperically [m/s]
rho_airport = 1.225;							# air density at runway altitude [kg/m^3]
sigma = rho_airport/rho0;						# Density ratio [-]
N_engines = 2;									# Number of engines [-]
c_v = 0.01199;									# Climb gradient divied by velocity [-]
Cl_max = array([2.9]);				            # Cl max values for assessment [-]
Cl_max_takeoff = array([2.5]);		            # Cl_max for take off [-]
A = array([9.5]);							# Aspect ratio [-]
W_S_max = 7000;									# Max Wing Loading value, change this value if you want to change the range of wing loading values you want to assess [N/m^2]
n = 2.5;										# Maximum load factor [-]
#%% ------------------------ Main ------------------------
T, p, rho_cruise, a = ISA_trop(h);
sigma_cruise = rho_cruise/rho0;			# Density ratio [-]
v_cruise = M*a;							# Cruise speed
PayloadRangeDiagram_JET(MTOW, OEW, Payload_max, 0.1, Max_Fuel_cap, (g, M, a, cj, L_D));

fig = plt.figure(figsize = (10, 8));
_ = W_S_stall(v_stall_landing, Cl_max, MLW/MTOW, fig);
_, _ = W_S_takeoff(Cl_max_takeoff, k, sigma, W_S_max, fig);
_, _ = W_S_cruise(A, Cd0, rho_cruise, sigma_cruise, v_cruise, W_S_max, fig);
_ = W_S_climb_grad(c_v, Cd0, A, e, W_S_max, N_engines, fig);
_, _ = W_S_maneuvering(n, Cd0, rho0, v_stall_landing, A, e, W_S_max, fig);

plt.grid(True);
plt.axis([0, W_S_max, 0, 1]);
plt.legend();
plt.xlabel("Wing Loading [N/m^2]");
plt.ylabel("Thrust Loading [-]");
plt.show();


from numpy import *
import matplotlib.pyplot as plt
from PayloadRangeDiagram import*
""" This program generates the wing and thrust loading diagram for a jet aircraft. All the variables in the
input block are the ones we can vary for the assessment. Please be careful with units as some emperical equations
only work with non SI units.
"""
#%% ---------------------- Constants ----------------------
g = 9.81;								# Gravitaional acceleration [m/s^2]
rho0 = 1.225;							# Sea level air Density [kg/m^3]
M = 0.78;								# Cruise mach number [-]
## ------------------------ INPUTS ------------------------
h = 11000;								# Cruise altitude [m]
MTOW = 80128*g;					# Max Take off Weight [N]
MLW = 0.86*MTOW                         # Max Landing Weight [N]
OEW = 44537*g;							# Operational empty weight [N]
Payload_max = 20*1000*g;						# Payload weight [N]
Max_Fuel_cap = 25*1000*g;
cj = 1.69895e-05;								# Specific fuel consumption [kg/Ns]
L_D = 16.7;										# Lift to drag ratio [-]
e = 0.794;								# Ozwald efficiency [-]
landdis = 1800;							# Landing distance [m]
Cd0 = 0.02139;							# Zero lift drag coefficient [-]
k = 565;								# Take off parameter from Raymer [retard units]
n_prop = 0.6                            # Propulsive efficiency (0.4-0.9)
v_stall_landing = sqrt(landdis/0.5847);			# Stall speed calcuated emperically [m/s]
rho_airport = 1.225;					# air density at runway altitude [kg/m^3]
sigma = rho_airport/rho0;				# Density ratio [-]
N_engines = 4;							# Number of engines [-]
Cl_max = array([2.2, 2.4, 2.6]);				# Cl max values for assessment [-]
Cl_max_takeoff = array([1.7, 1.9, 2.1]);	# Cl_max for take off [-]
A = array([9, 10, 11]);					# Aspect ratio [-]
W_S_max = 7000;							# Max Wing Loading value, change this value if you want to change the range of wing loading values you want to assess [N/m^2]
n = 2.5;								# Maximum load factor [-]
#%% ---------------------- Functions ----------------------

def W_S_stall(v_stall_landing, Cl_max_landing, MLW_MTOW, fig):
	""" Function to compute the wing loading
		for stall speeds
	Input:
		v_stall_landing = Stall speed in landing configuration (float) [m/s]
		Cl_max_landing = Max Cl value in landing configuration (float or array) [-]
		MLW_MTOW = Max landing weight divied by MTOW [-]
		fig = Figure handel to plot to (handel)
	Output:
		W_load_Stall = Wing Loading (float or array) [N/m^2]
	"""
	W_load_Stall = 0.5*1.225*(v_stall_landing)**2*Cl_max_landing/(MLW_MTOW);
	plt.figure(fig.number);
	try:									# So that function can except both float and array values for Cl
		for i in range(len(W_load_Stall)):
			plt.plot([W_load_Stall[i], W_load_Stall[i]], [0, 1], label = "Cl land = " + str(Cl_max_landing[i]));
	except:
		plt.plot([W_load_Stall, W_load_Stall], [0, 1], label = "Cl land = " + str(Cl_max_landing));
	return W_load_Stall;

def W_S_takeoff(Cl_max_takeoff, k, sigma, W_S_max, fig):
	""" Function assess take off performance
	Ref Slide 26 Lecture 3 from ADSEE 1. Also refer to "TakeoffPerformance.png"
	for the value of k.
	Input:
		Cl_max_takeoff = Max Cl value at take-off [-] (float or array)
		k = Take off parameter, directly read from "TakeoffPerformance.png" [retard units but already acounted for in this function]
		sigma = Density ratio [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to (handel)
	Output:
		x = Range of wing loading values [N/m^2] (array)
		y = Range of thrust loading values [N/N] (array)
	"""
	k = k * 4.45**2 / (745.6*0.3048**2);										# Convert from BHP/lbs to Watt/N (as a sane person whould do)
	coeff = (Cl_max_takeoff/1.21) *k*sigma;
	x = linspace(1, W_S_max, 100);						# Power Loading, W/P [N/W]
	# plot to fig
	plt.figure(fig.number);
	try:												# So that function can except both float and array values for Cl
		y = coeff.reshape(len(coeff), 1) / x;				# Wing Loading, W/S [N/m^2]
		for i in range(len(coeff)):
			plt.plot(x, y[i, :], label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(round_(Cl_max_takeoff[i], 2)));
	except:
		y = x*coeff;									# Wing Loading, W/S [N/m^2]
		plt.plot(x, y, label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(round_(Cl_max_takeoff, 2)));
	return x, y;

def W_S_cruise(A, CD0, rho_cruise, sigma, v_cruise, n_prop, W_S_max, fig):
	""" Function to assess Wing and Thrust loading during cruise
	Input:
		A = Aspect ratio [-] (float or array)
		CD0 = CD0 zero lift drag coefficient [-] (float)
		rho_cruise = air density at cruise altitude [kg/m^3] (float)
		v_cruise = cruise speed [m/s] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to (handel)
	Output:
		x = Range of wing loading values [N/m^2] (array)
		y = Range of thrust loading values [N/N] (array)
	"""
	throttle = 0.9;
	cru_weight = 0.8;
	rel_rho = 1/sigma;
	x = linspace(1, W_S_max, 1000);
	plt.figure(fig.number);
	try:
		y = throttle/cru_weight*n_prop*rel_rho**(3/4)*((CD0*0.5*rho_cruise*v_cruise**3/(cru_weight*x)) +
										((cru_weight*x)*1/(pi*A.reshape(len(A), 1)*e*0.5*rho_cruise*v_cruise)))**(-1);
		for i in range(len(y)):
			plt.plot(x, y[i, :], "--", label = "Aspect Ratio for cruise = " + str(A[i]));
	except:
		y = throttle/cru_weight*n_prop*rel_rho**(3/4)*((CD0*0.5*rho_cruise*v_cruise**3/(cru_weight*x)) +
										((cru_weight*x)*1/(pi*A*e*0.5*rho_cruise*v_cruise)))**(-1);
		plt.plot(x, y, "--", label = "Aspect Ratio for cruise = " + str(A));
	return x, y;

def W_S_climb_grad(Cl_max, Cl_max_takeoff, n_prop, Cd0, A, e, W_S_max, N_engines, fig):
	""" This function computes the Wing Loading to meet the set climb gradient
		requirment, c_v in a one engine inoerative case
	Input:
		Cl_max = maximum lift coefficient at landing
		Cl_max_takeoff = max lift coefficient at take-off
		n_prop = propulsive efficiency, between 0.4 and 0.9 heavily dependent on velocity (P_a = n_prop * P_br)
		Cd0 = zero lift drag coefficient [-] (float)
		A = Aspect ratio [-] (float or array)
		e = Ozwald efficiency [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		N_engines = Number of engines
		fig = Figure handel to plot to (handel)
	Output:
		y = Range of power loading values [N/W] (array)
	"""
	Cd0_goaround = Cd0 + 0.065 + 0.025;			 #source ADSEE 13 slide 13 (Roskam)
	e_goaround = e + 0.1;
	CL_goaround = Cl_max/1.21
	CL_goaround = CL_goaround.reshape(len(CL_goaround), 1)
	CD_goaround = Cd0_goaround + CL_goaround**2/(pi*A.reshape(len(A), 1)*e_goaround)
	CD_goaround = CD_goaround.reshape(len(CD_goaround), 1)

	Cd0_OEI = Cd0 + 0.015;
	e_OEI = e + 0.05;
	CL_OEI = Cl_max_takeoff/1.21
	CL_OEI = CL_OEI.reshape(len(CL_OEI), 1)
	CD_OEI = Cd0_OEI + CL_OEI**2/(pi*A.reshape(len(A), 1)*e_OEI)
	CL_OEI = CD_OEI.reshape(len(CD_OEI), 1)

	if N_engines == 2:
		c_v_OEI = 0.024; #CS25.121
		c_v_goaround = 0.032 #CS25.121
	elif N_engines == 3:
		c_v_OEI = 0.027; #CS25.121
		c_v_goaround = 0.032 #CS25.121
	elif N_engines > 3:
		c_v_OEI = 0.03 #CS25.121
		c_v_goaround = 0.032 #CS25.121
	else:
		print("wrong definition of number of engines")

	x = linspace(1, W_S_max, 1000)
	plt.figure(fig.number);




	try:
		y_goaround = n_prop/(sqrt(x)*(c_v_goaround+CD_goaround/CL_goaround)*sqrt(2/(1.225*CL_goaround)))
		y_OEI = ((N_engines-1)/N_engines)*n_prop/(sqrt(x)*(c_v_OEI+CD_OEI/CL_OEI)*sqrt(2/(1.225*CL_OEI)))

		for i in range(len(y_goaround)):
			plt.plot(x, y_goaround[i, :], "--", label = "CL_max_land = " + str(Cl_max[i]));
			plt.plot(x, y_OEI[i, :], "--", label = "CL_max_takeoff = " + str(Cl_max_takeoff[i]));
	except:
		y_goaround = n_prop/(sqrt(x)*(c_v_goaround+CD_goaround/CL_goaround)*sqrt(2/(1.225*CL_goaround)))
		y_OEI = ((N_engines-1)/N_engines)*n_prop/(sqrt(x)*(c_v_OEI+CD_OEI/CL_OEI)*sqrt(2/(1.225*CL_OEI)))

		plt.plot(x, y_goaround[i, :], "--", label = "CL_max_land = " + str(Cl_max[i]));
		plt.plot(x, y_OEI[i, :], "--", label = "CL_max_takeoff = " + str(Cl_max_takeoff[i]));
	return y_goaround;

def W_S_maneuvering(n, Cd0, rho, v_stall_landing, A, e, W_S_max, fig):
	""" This function computes the Wing Loading to meet the maneuvering
		requirment, given a max load factor of n
	Input:
		n = Max load factor [-] (float)
		Cd0 = zero lift drag coefficient [-] (float)
		v = Velocity [m/s] (float)
		A = Aspect ratio [-] (float or array)
		e = Ozwald efficiency [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to (handel)
	Output:
		y = Range of thrust loading values [N/N] (array)
	"""
	v_stall_clean = 1.2 * v_stall_landing		# not correct, awaiting advice from course lecturer
	v = v_stall_clean*sqrt(n);					# this is the manouvre speed
	x = linspace(1, W_S_max, 1000);
	plt.figure(fig.number);
	try:
		y = Cd0*0.5*rho*v**2/x + x*n**2/(pi*A.reshape(len(A), 1)*e*0.5*rho*v**2);
		for i in range(len(y)):
			plt.plot(x, y[i], label = "Aspect Ratio for maneuvering = " + str(A[i]));
	except:
		y = Cd0*0.5*rho*v**2/x + x*n**2/(pi*A*e*0.5*rho*v**2);
		plt.plot(x, y, label = "Aspect Ratio for maneuvering = " + str(A));
	return x, y;
#%% ---------------------- Main ----------------------
T, p, rho_cruise, a = ISA_trop(h);
sigma_cruise = rho_cruise/rho0;			# Density ratio [-]
v_cruise = M*a;							# Cruise speed [m/s]
PayloadRangeDiagram_JET(MTOW, OEW, Payload_max, 0.1, Max_Fuel_cap, (g, M, a, cj, L_D));



fig = plt.figure(figsize = (10, 8));
MLW_MTOW = MLW / MTOW
_ = W_S_stall(v_stall_landing, Cl_max, MLW_MTOW, fig);
_, _ = W_S_takeoff(Cl_max_takeoff, k, sigma, W_S_max, fig);
_, _ = W_S_cruise(A, Cd0, rho_cruise, sigma_cruise, v_cruise, n_prop, W_S_max, fig);
_ = W_S_climb_grad(Cl_max, Cl_max_takeoff, n_prop, Cd0, A, e, W_S_max, N_engines, fig);
#_, _ = W_S_maneuvering(n, Cd0, rho0, v_stall_landing, A, e, W_S_max, fig);

plt.grid(True);
plt.axis([0, W_S_max, 0, 1]);
plt.legend();
plt.xlabel("Wing Loading [N/m^2]");
plt.ylabel("Power Loading [N/W]");
plt.show();


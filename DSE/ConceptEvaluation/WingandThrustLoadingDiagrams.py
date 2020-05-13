from numpy import *
import matplotlib.pyplot as plt
""" This program generates the wing and thrust loading diagram for a jet aircraft. All the variables in the 
input block are the ones we can vary for the assessment. Please be careful with units as some emperical equations
only work with non SI units. 
"""
#%% ---------------------- Constants ----------------------
g = 9.81;								# Gravitaional acceleration [m/s^2]
rho0 = 1.225;							# Sea level air Density [kg/m^3]
M = 0.78;								# Cruise mach number [-]
# ------------------------ INPUTS ------------------------
h = 10000;								# Cruise altitude [m]
MTOW = 77530.96818*g;					# Max Take off Weight [N]
e = 0.794;								# Ozwald efficiency [-]
landdis = 1600;							# Landing distance [m]
Cd0 = 0.0156;							# Zero lift drag coefficient [-]
k = 120;								# Take off parameter from "TakeOffPerformance.png" [retard units]
v_stall = sqrt(landdis/0.5847);			# Stall speed calcuated emperically [m/s]
rho_airport = 1.225;					# air density at runway altitude [kg/m^3]
sigma = rho_airport/rho0;				# Density ratio [-]
N_engines = 2;							# Number of engines [-]
c_v = 0.023993;							# Climb gradient divied by velocity [-]
Cl_max = array([1.7, 2.5]);				# Cl max values for assessment [-]
Cl_max_take_off = array([1.5, 1.7]);	# Cl_max for take off [-]
A = array([6, 8.4, 9.5]);				# Aspect ratio [-]
W_S_max = 5000;							# Max Wing Loading value, change this value if you want to change the range of wing loading values you want to assess [N/m^2]
n = 2.5;								# Maximum load factor [-]
MLW_MTOW = 0.86;						# Max landing weight divided by max take off weight [-]
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
	a = sqrt(1.4*287*T);
	return T, p, rho, a;

def W_S_stall(v_land, rho, Cl_max, MLW_MTOW, fig):
	""" Function to compute the wing loading
		for stall speeds
	Input:
		v_land = 1.1*Stall speed (float) [m/s]
		rho = air density (float) [kg/m^3]
		Cl_max = Max Cl value (float or array) [-]
		MLW_MTOW = Max landing weight divided by max take off weight [-]
		fig = Figure handel to plot to (handel)
	Output:
		W_load_Stall = Wing Loading (float or array) [N/m^2]
	"""
	W_load_Stall = 0.5*rho*v_land**2*Cl_max/MLW_MTOW;
	plt.figure(fig.number);
	try:									# So that function can except both float and array values for Cl
		for i in range(len(W_load_Stall)):
			plt.plot([W_load_Stall[i], W_load_Stall[i]], [0, 1], label = "Cl land = " + str(Cl_max[i]));
	except:
		plt.plot([W_load_Stall, W_load_Stall], [0, 1], label = "Cl land = " + str(Cl_max));
	return W_load_Stall;

def W_S_takeoff(Cl_max, k, sigma, W_S_max, fig):
	""" Function assess take off performance
	Ref Slide 26 Lecture 3 from ADSEE 1. Also refer to "TakeoffPerformance.png"
	for the value of k.
	Input:
		Cl_max = Max Cl value [-] (float or array)
		k = Take off parameter, directly read from "TakeoffPerformance.png" [retard units but already acounted for in this function]
		sigma = Density ratio [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to (handel)
	Output:
		x = Range of wing loading values [N/m^2] (array)
		y = Range of thrust loading values [N/N] (array)
	"""
	k *= 47.880172;										# Convert from retard units to SI units (as a sane person whould do)
	coeff = 1/(Cl_max*k*sigma);
	x = linspace(1, W_S_max, 100);						# Thrust Loading, T/W [N/N]
	# plot to fig
	plt.figure(fig.number);
	try:												# So that function can except both float and array values for Cl
		y = x*coeff.reshape(len(coeff), 1);				# Wing Loading, W/S [N/m^2]
		for i in range(len(coeff)):
			plt.plot(x, y[i, :], label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(round_(Cl_max[i], 2)));
	except:
		y = x*coeff;									# Wing Loading, W/S [N/m^2]
		plt.plot(x, y, label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(round_(Cl_max, 2)));
	return x, y;

def W_S_cruise(A, CD0, rho_cruise, sigma, v_cruise, W_S_max, fig):
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
		y = cru_weight/throttle*rel_rho**(3/4)*((CD0*0.5*rho_cruise*v_cruise**2/(cru_weight*x)) +
										((cru_weight*x)*1/(pi*A.reshape(len(A), 1)*e*0.5*rho_cruise*v_cruise**2)));
		for i in range(len(y)):
			plt.plot(x, y[i, :], "--", label = "Aspect Ratio for cruise = " + str(A[i]));
	except:
		y = cru_weight/throttle*rel_rho**(3/4)*((CD0*0.5*rho_cruise*v_cruise**2/(cru_weight*x)) +
										((cru_weight*x)*1/(pi*A*e*0.5*rho_cruise*v_cruise**2)));
		plt.plot(x, y, "--", label = "Aspect Ratio for cruise = " + str(A));
	return x, y;

def W_S_climb_grad(c_v, Cd0, A, e, W_S_max, N_engines, fig):
	""" This function computes the Wing Loading to meet the set climb gradient
		requirment, c_v in a one engine inoerative case
	Input:
		c_v = required climb gradent divied by velocity [rad*s/m] (float)
		Cd0 = zero lift drag coefficient [-] (float)
		A = Aspect ratio [-] (float or array)
		e = Ozwald efficiency [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		N_engines = Number of engines
		fig = Figure handel to plot to (handel)
	Output:
		y = Range of thrust loading values [N/N] (array)
	"""
	y = (c_v + 2*sqrt(Cd0/(pi*A*e)));
	y_OEI = N_engines/(N_engines - 1)*(c_v + 2*sqrt(Cd0/(pi*A*e)));
	plt.figure(fig.number);
	try:
		for i in range(len(y)):
			plt.plot([0, W_S_max], [y[i], y[i]], label = "Aspect Ratio climb gradient = " + str(A[i]));
			plt.plot([0, W_S_max], [y_OEI[i], y_OEI[i]], label = "Aspect Ratio (one engine inoperative) = " + str(A[i]));
	except:
		plt.plot([0, W_S_max], [y, y], label = "Aspect Ratio climb gradient = " + str(A));
		plt.plot([0, W_S_max], [y_OEI, y_OEI], label = "Aspect Ratio (one engine inoperative) = " + str(A));
	return y;

def W_S_maneuvering(n, Cd0, rho, v, A, e, W_S_max, fig):
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
_, _, rho_cruise, a = ISA_trop(h);
sigma_cruise = rho_cruise/rho0;			# Density ratio [-]
v_cruise = M*a;							# Cruise speed [m/s]
v_land = v_stall*1.1;					# landing speed [m/s]
fig = plt.figure(figsize = (10, 8));
_ = W_S_stall(v_land, rho_airport, Cl_max, MLW_MTOW, fig);
_, _ = W_S_takeoff(Cl_max_take_off, k, sigma, W_S_max, fig);
_, _ = W_S_cruise(A, Cd0, rho_cruise, sigma_cruise, v_cruise, W_S_max, fig);
_ = W_S_climb_grad(c_v, Cd0, A, e, W_S_max, N_engines, fig);
_, _ = W_S_maneuvering(n, Cd0, rho_cruise, v_cruise, A, e, W_S_max, fig);

plt.grid(True);
plt.axis([0, W_S_max, 0, 1]);
plt.legend();
plt.xlabel("Wing Loading [N/m^2]");
plt.ylabel("Thrust Loading [-]");
plt.show();
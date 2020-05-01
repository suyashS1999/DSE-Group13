from numpy import *
import matplotlib.pyplot as plt

#%% ---------------------- Constants ----------------------
g = 9.81;						# [m/s^2]
rho0 = 1.225;					# Sea level air Density [kg/m^3]
# ------------------------ INPUTS ------------------------
MTOW = 77530.96818;				# [kg]
S = 257.;						# Wing Area [m^2]
e = 0.794;						# Ozwald efficiency [-]
W_S_takeoffdis = 2100;			# Take off distance [ft]
landdis = 1600;					# Landing distance [ft]
Cd0 = 0.0156;					# Cd0 value [-]
k = 120;						# Take off parameter from "W_S_takeoffPerformance.png" [retard units]
v_stall = sqrt(landdis/0.5847);	# Stall speed calcuated emperically [m/s]
rho_airport = 1.225;			# air density at runway altitude [kg/m^3]
rho_cruise = 0.33;				# air density at cruise altitude [kg/m^3]
sigma = rho_airport/rho0;		# Density ratio [-]
sigma_cruise = rho_cruise/rho0	# Density ratio [-]
v_cruise = 232.6;				# Cruise speed [m/s]
c_v = 0.023993;					# Climb gradient divied by velocity [s/m]
Cl_max = array([1.7, 2.5]);		# Cl max values for assessment [-]
A = array([6, 8.4, 9.5]);		# Aspect ratio [-]
W_S_max = 5000;					# Max Thrust Loading value [N/m^2]
#%% ---------------------- Functions ----------------------
def W_S_stall(v_stall, rho, Cl_max, fig):
	""" Function to compute the wing loading
		for stall speeds
	Input:
		v_stall = Stall speed (float) [m/s]
		rho = air density (float) [kg/m^3]
		Cl_max = Max Cl value (float) [-]
	Output:
		W_load_Stall = Wing Loading (float) [N/m^2]
	"""
	W_load_Stall = 0.5*rho*v_stall**2*Cl_max;
	plt.figure(fig.number);
	try:									# So that function can except both float and array values for Cl
		for i in range(len(W_load_Stall)):
			plt.plot([W_load_Stall[i], W_load_Stall[i]], [0, 1], label = "Cl land = "+ str(Cl_max[i]));
	except:
		plt.plot([W_load_Stall, W_load_Stall], [0, 1], label = "Cl land = "+ str(Cl_max));
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
		fig = Figure handel to plot to
	Output:
		x = Range of thrust loading values [N/N] (array)
		y = Range of wing loading values for the given Take Off Performance [N/m^2] (array)
	"""
	k *= 47.880172;										# Convert from retard units to SI units (as a sane person whould do XD)
	coeff = 1/((Cl_max/1.21)*k*sigma);
	x = linspace(1, W_S_max, 100);						# Thrust Loading, T/W [N/N]
	# plot to fig
	plt.figure(fig.number);
	try:												# So that function can except both float and array values for Cl
		y = x*coeff.reshape(len(coeff), 1);				# Wing Loading, W/S [N/m^2]
		for i in range(len(coeff)):
			plt.plot(x, y[i, :], label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(Cl_max[i]));
	except:
		y = x*coeff;									# Wing Loading, W/S [N/m^2]
		plt.plot(x, y, label = "TOP = " + str(round_(k, 5)) + " , CL take off = " + str(Cl_max));
	return x, y;

def W_S_cruise(A, CD0, rho_cruise, sigma, v_cruise, W_S_max, fig):
	""" Function to assess Wing and Thrust loading during cruise
	Input:
		A = Aspect ratio [-] (float or array)
		CD0 = CD0 zero lift drag coefficient [-] (float)
		rho_cruise = air density at cruise altitude [kg/m^3] (float)
		v_cruise = cruise speed [m/s] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to
	Output:
		x = Range of thrust loading values [N/N] (array)
		y = Range of wing loading values for Cruise Performance [N/m^2] (array)
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

	return x, y

def W_S_climb_grad(c_v, Cd0, A, e, W_S_max, fig):
	""" This function computes the Wing Loading to meet the set climb gradient
		requirment, c_v in a one engine inoerative case
	Input:
		c_v = required climb gradent divied by velocity [rad*s/m] (float)
		Cd0 = zero lift drag coefficient [-] (float)
		A = Aspect ratio [-] (float or array)
		e = Ozwald efficiency [-] (float)
		W_S_max = maximum wing loading to be expected [N/m^2] (float)
		fig = Figure handel to plot to
	"""
#Design for one engine inoperative
	y = (c_v + 2*sqrt(Cd0/(pi*A*e)));
	plt.figure(fig.number);
	try:
		for i in range(len(y)):
			plt.plot([0, W_S_max], [y[i], y[i]], label = "Aspect Ratio (one engine inoperative) = "+ str(A[i]));
	except:
		plt.plot([0, W_S_max], [y, y], label = "Aspect Ratio (one engine inoperative) = "+ str(A));

#%% ---------------------- Main ----------------------
fig = plt.figure(figsize = (10, 8));
_ = W_S_stall(v_stall, rho_airport, Cl_max, fig);
_, _ = W_S_takeoff(Cl_max, k, sigma, W_S_max, fig);
_, _ = W_S_cruise(A, Cd0, rho_cruise, sigma_cruise, v_cruise, W_S_max, fig);
_ = W_S_climb_grad(c_v, Cd0, A, e, W_S_max, fig);

plt.grid(True);
plt.axis([0, 5000, 0, 1]);
plt.legend();
plt.xlabel("Wing Loading [W/S]");
plt.ylabel("Thrust Loading [T/W]");
plt.show();
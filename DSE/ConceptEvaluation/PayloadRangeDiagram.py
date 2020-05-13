# -*- coding: utf-8 -*-
"""
Used in the conceptual stage to estimate main parameters for concept selection
Outputs are
Payload range diagram
Weight
Volume
Cost 
"""

from numpy import*
from matplotlib import pyplot as plt

# ---------------------- Constants ----------------------
g = 9.81;						# [m/s^2]
#M = 0.78;						# [-]
##%% ------------------------ INPUTS ------------------------
#MTOW = 76003.00497;				# [kg]
#OEW = 41195.17484;				# [kg]
#Payload_max = 20*1000;			# [kg]
#Max_Fuel_cap = 25*1000;			# [kg]
## ---------------------- Variables ----------------------
#h = 10000;						# [m]
#cj = 1.69895e-05;				# [kg/Ns]
#L_D = 20;						# [-]
#Weight_frac = 1/0.857812138;			# [-]

#%% ------------------------ Functions ------------------------
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

def RangeJet(g, M, a, cj, L_D, Weight_frac):
	""" This function computes the range that can be obtained by a 
		jet aircraft
	Input:
		V = cruise velocity [m/s]
		cj = Thrust specific fuel consumption [kg/Ns]
		L_D = Lift to Drag ratio [-]
		Weight_frac = weight fraction before and after cruise phase [-]
	Output:
		R = Range [m]
	"""
	V = M*a;
	R = (V/(cj*g))*L_D*log(Weight_frac);
	return R;

def PayloadRangeDiagram_JET(MTOW, OEW, Payload_max, reserve_fuel_frac, fuel_ascent_descent_frac, Max_Fuel_cap, Range_func_args):
	""" This function generates the payload range diagram 
		for a given JET aircraft.
	Input:
		MTOW = Max Take Off Weight [kg or N]
		OEW = Operational Empty Weight [kg or N]
		Payload_max = Maximum Payload [kg or N]
		reserve_fuel_frac = Fraction of the total fuel that is reserve fuel [-]
		fuel_ascent_descent_frac = Fuel fraction used during ascending and descending [-]
		Max_Fuel_cap = Maximum fuel capacity [kg]
		Range_func_args = Arguments to pass onto Range function
	Output:
		Payload Range Diagram Plot
	"""
	fig = plt.figure(figsize = (14, 8));
	# Harmonic Range
	MZFW = OEW + Payload_max;
	Fuel_Weight = MTOW - MZFW;	Reserve_fuel = reserve_fuel_frac*Fuel_Weight;	Trip_fuel = Fuel_Weight - Reserve_fuel;
	Fuel_ac_dc = fuel_ascent_descent_frac*Fuel_Weight;
	MZFW_inc_res_fuel = MZFW + Reserve_fuel;
	MZFW_inc_res_ac_dc_fuel = MZFW + Reserve_fuel + Fuel_ac_dc;
	fuel_arr = linspace(0, MTOW - MZFW_inc_res_ac_dc_fuel, 100);
	R_Harmonic = RangeJet(*Range_func_args, (MZFW_inc_res_ac_dc_fuel + fuel_arr)/MZFW_inc_res_ac_dc_fuel)/1000;
	
	plt.plot(R_Harmonic, MZFW*ones_like(R_Harmonic));
	plt.plot(R_Harmonic, MZFW_inc_res_fuel*ones_like(R_Harmonic), ":");
	plt.plot(R_Harmonic, MZFW_inc_res_ac_dc_fuel*ones_like(R_Harmonic), ":");
	plt.annotate("MZFW", xy = (R_Harmonic[int(len(R_Harmonic)/2)], MZFW), xytext = (R_Harmonic[int(len(R_Harmonic)/2)], MZFW*(0.9)), 
			  ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	plt.plot(R_Harmonic, MZFW_inc_res_ac_dc_fuel + fuel_arr);
	plt.annotate("Reserve Fuel", xy = (R_Harmonic[int(len(R_Harmonic)/2) + 20], MZFW_inc_res_fuel), 
			xytext = (R_Harmonic[int(len(R_Harmonic)/2) + 20], MZFW*(0.99)), 
			ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "<->"));
	plt.plot([R_Harmonic[-1], R_Harmonic[-1]], [OEW, MTOW], color = "r", label = "Harmonic Range");
	plt.fill_between(R_Harmonic, MZFW*ones_like(R_Harmonic), OEW, facecolor = "orange", color = "blue", alpha = 0.2, label = "Payload");

	# Design Range
	payload_arr = linspace(MZFW_inc_res_ac_dc_fuel, (MTOW - Max_Fuel_cap + Reserve_fuel + Fuel_ac_dc), 100);
	R_Design = RangeJet(*Range_func_args, MTOW/payload_arr)/1000;

	plt.plot(R_Design, payload_arr, ":");
	plt.plot(R_Design, MTOW*ones_like(R_Design));
	plt.plot(R_Design, payload_arr - Reserve_fuel - Fuel_ac_dc);
	plt.plot(R_Design, payload_arr - Fuel_ac_dc, ":");
	idx = int(len(R_Design)/2) - 15;
	plt.annotate("Trip Fuel", xy = (R_Design[idx], payload_arr[idx] - Fuel_ac_dc), 
		xytext = (R_Design[idx], (MTOW - payload_arr[idx])/2 + payload_arr[idx]), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	plt.annotate("", xy = (R_Design[idx], MTOW), 
		xytext = (R_Design[idx], (MTOW - payload_arr[idx])/2 + payload_arr[idx]*(1.01)), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	idx += 20;
	plt.annotate("Fuel ascent/ \n descent", xy = (R_Design[idx], payload_arr[idx] - Fuel_ac_dc), 
		xytext = (R_Design[idx], payload_arr[idx]), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	plt.annotate("", xy = (R_Design[idx], payload_arr[idx]), 
		xytext = (R_Design[idx], payload_arr[idx]*(0.99)), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	plt.plot([R_Design[-1], R_Design[-1]], [OEW, MTOW], color = "c", label = "Max Design Range");
	plt.fill_between(R_Design, payload_arr - Reserve_fuel - Fuel_ac_dc, OEW, facecolor = "orange", color = "blue", alpha = 0.2);

	# Ferry Range
	payload_ferry_arr = linspace(MTOW - Max_Fuel_cap + Reserve_fuel + Fuel_ac_dc, OEW + Reserve_fuel + Fuel_ac_dc, 100);
	TOW = (MTOW - (payload_ferry_arr[0] - payload_ferry_arr));
	R_Ferry = RangeJet(*Range_func_args, TOW/(payload_ferry_arr))/1000;

	plt.plot(R_Ferry, TOW);
	plt.plot(R_Ferry, payload_ferry_arr, ":");
	plt.plot(R_Ferry, payload_ferry_arr - Reserve_fuel - Fuel_ac_dc);
	plt.plot(R_Ferry, payload_ferry_arr - Fuel_ac_dc, ":");
	plt.fill_between(R_Ferry, payload_ferry_arr - Reserve_fuel - Fuel_ac_dc, OEW, facecolor = "orange", color = "blue", alpha = 0.2);
	idx = int(len(R_Ferry)/2);
	plt.annotate("Maximum Fuel Capacity", xy = (R_Ferry[idx], payload_ferry_arr[idx] - Reserve_fuel - Fuel_ac_dc), 
		xytext = (R_Ferry[idx], (TOW[idx] - payload_ferry_arr[idx])/2 + payload_ferry_arr[idx]), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));
	plt.annotate("", xy = (R_Ferry[idx], TOW[idx]), 
		xytext = (R_Ferry[idx], (TOW[idx] - payload_ferry_arr[idx])/2 + payload_ferry_arr[idx]*(1.01)), 
		ha = "center", arrowprops = dict(facecolor = "black", arrowstyle = "->"));

	R = concatenate((R_Harmonic, R_Design, R_Ferry));
	plt.plot(R, MTOW*ones_like(R), '--', label = "MTOW");
	plt.plot(R, OEW*ones_like(R), '--', label = "OEW");
	plt.plot([R[-1], R[-1]], [TOW[-1], OEW], "k-", label = "Ferry Range");

	plt.xlabel("Range [km]");
	plt.ylabel("Weight [N]");
	plt.grid(True);
	plt.legend(loc = 1);
	plt.tight_layout();
	plt.xlim(xmin = 0);
	#plt.show();
	return 0;

#%% ------------------------ MAIN -----------------------
#T, p, rho, a = ISA_trop(h);
#R = RangeJet(g, M, a, cj, L_D, Weight_frac);
#PayloadRangeDiagram_JET(MTOW, OEW, Payload_max, 0.1, Max_Fuel_cap, (g, M, a, cj, L_D));

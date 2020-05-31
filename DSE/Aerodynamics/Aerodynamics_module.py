from numpy import*
from matplotlib import pyplot as plt
""" Code for airfoil analysis. Takes input from XFLR5 data and corrects for compressible
flow and predicts Mach number over airfoil. Note, the code used linearized theory for 
compressible flow thus this analysis is only valid for small pertubations ie. small AOA.
"""
def XFLR5PolarDataextraction(path_polars, path_X):
	""" This function extracts the simulation data from 
		the .txt files generated with XFLR5 software for 2D airfoil simulation
	Input:
		path_polars = string with a "r" in the front giving the path to the location of the txt file containing the polars
		path_X = string with a "r" in the front giving the path to the location of the txt file containing the chord locations
	Output:
		polars = array containing the AOA, Cl, Cd, Cm and Cp values for each AOA (array)
		X = chord location of the Cp values (array)
		Cp_start_idx = array containing the indicies of where the first Cp value is for each AOA (array)
	"""
	polars = array(open(path_polars).read().splitlines());
	unwantedtxt = [polars[0], polars[1], polars[2], polars[3], polars[4], polars[6], polars[7]];
	X = genfromtxt(path_X)[:, 0];

	for i in range(len(unwantedtxt)):
		idx = where(polars == unwantedtxt[i]);
		if i == len(unwantedtxt) - 1:
			Cp_start_idx = idx[0] - range(0, len(idx[0]));
		polars = delete(polars, idx, 0);

	polars = ([x.split() for x in polars]);
	polars = array([[float(j) for j in i] for i in (polars)]);
	return polars, X, Cp_start_idx;

def ApplyCompressibilityCorrection(M, C):
	""" This function applies the Prandtl Glauert correction for 
		a given free-stram mach number to the quantity C
	Input:
		M = Freestram Mach number [-] (float or int)
		C = any aerodynamic quantity that needs compressibility correction (array or float or int)
	Output:
		C = corrected quantity (same a input C)
	"""
	beta = (sqrt(1 - M**2));
	coeff = 1/(beta)# + 0.5*(1 - beta)*C);	# + (M**2/(1 + beta))*C/2);
	return C*coeff;

def ComputeLocalPressure(M, Cp, p_inf, gamma, ratio):
	""" This function computes the value of the static pressure or the 
		pressure ratio given the Cp value at the point
	Input:
		M = FreeStream Mach number [-] (float or int)
		Cp = Coefficient of pressure [-] (array or float)
		p_inf = Freestream static pressure [Pa] (float)
		gamma = ratio of specific heats [-] (float)
		ratio = Boolien, if ratio == True: return Pressure ratio if False: return static Pressure
	Output:
		P = Pressure value or pressure ratio
	"""
	P = ((Cp*gamma*M**2)/2 + 1)*p_inf;
	P_tot =	amax(P);				#p_inf*(1 + (gamma - 1)/2*M**2)**(gamma/(gamma - 1));
	if ratio == True:
		return P_tot/P;
	elif ratio == False:
		return P;

def IsentropicRelation(gamma, ratio_value, P_or_T_ratio):
	""" This function computes the value of the Mach number 
		using Isentropic relations
	Input:
		gamma = ratio of specific heats [-] (float)
		ratio_value = Pressure ratio or Temperature ratio [-] (float or array)
		P_or_T_ratio = Specify whether ratio_values is Pressure ratio or Temperature ratio [-] (string)
	Output:
		M = Mach number
	"""
	if P_or_T_ratio == "Pressure":
		ratio_value = sign(ratio_value)*absolute(ratio_value)**((gamma - 1)/gamma);
	M = sqrt((ratio_value - 1)*2/(gamma - 1));
	return M;

#%% Input data
# These paths are specific to my PC, change path to your own directory, and don't forget the "r" at the start
f_Polars = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\Polars2.txt";
f_X = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\x_col.txt";
M = 0.78;				# Freestream Mach 
p_inf = 100000.;		# Freestream staic pressure
gamma = 1.4;			# ratio of specific heats

aoa = 2;				# Specify the aoa you want to analyse
print("Specified Angle of attack: {} degrees".format(aoa));
polars, X, Cp_start_idx = XFLR5PolarDataextraction(f_Polars, f_X);
data = vstack((polars)[Cp_start_idx - 1]);
AOA = data[:, 0];
Cl = data[:, 2];		Cl_comp = ApplyCompressibilityCorrection(M, Cl);
Cd = data[:, 1];		Cd_comp = ApplyCompressibilityCorrection(M, Cd);
Cm = data[:, 3];		Cm_comp = ApplyCompressibilityCorrection(M, Cm);

CP_i = zeros((len(X), len(AOA)));
CP_v = zeros((len(X), len(AOA)));
for i in range(len(AOA)):
	idx = Cp_start_idx[i];
	CP_i[:, i] = vstack(polars[idx: idx + len(X)])[:, 0];
	CP_v[:, i] = vstack(polars[idx: idx + len(X)])[:, 1];

CP = CP_v;

#%% Plots
fig = plt.figure(figsize = (16, 16));
ax1 = plt.subplot(2, 2, 1);
ax1.plot(AOA, Cl, label = "Incompressible");
ax1.plot(AOA, Cl_comp, 'r', label = "Compressible M =" + str(M));
ax1.set_xlabel("alpha [degrees]");
ax1.set_ylabel("Cl [-]");

ax2 = plt.subplot(2, 2, 2);
ax2.plot(AOA, Cd, label = "Incompressible");
ax2.plot(AOA, Cd_comp, 'r', label = "Compressible M =" + str(M));
ax2.set_xlabel("alpha [degrees]");
ax2.set_ylabel("Cd [-]");

ax3 = plt.subplot(2, 2, 3);
ax3.plot(AOA, Cm, label = "Incompressible");
ax3.plot(AOA, Cm_comp, 'r', label = "Compressible M =" + str(M));
ax3.set_xlabel("alpha [degrees]");
ax3.set_ylabel("Cm [-]");

ax4 = plt.subplot(2, 2, 4);
ax4.plot(Cd, Cl, label = "Incompressible");
ax4.plot(Cd_comp, Cl_comp, 'r', label = "Compressible M =" + str(M));
ax4.set_xlabel("Cd [-]");
ax4.set_ylabel("Cl [-]");
ax1.grid(True);
ax2.grid(True);
ax3.grid(True);
ax4.grid(True);
ax1.legend();
ax2.legend();
ax3.legend();
ax4.legend();

req_idx = where(AOA == aoa);
cp = CP[:, req_idx].reshape(1, -1)[0];
cp_comp = ApplyCompressibilityCorrection(M, cp);
P_local = ComputeLocalPressure(M, cp_comp, p_inf, gamma, True);
M_local = IsentropicRelation(gamma, P_local, "Pressure");

fig = plt.figure(figsize = (16, 8));
ax1 = plt.subplot(1, 2, 1);
ax1.plot(X, cp, 'bx-', label = "Incompressible");
ax1.plot(X, cp_comp, 'rx-', label = "Compressible M =" + str(M));
ax1.invert_yaxis();
ax1.grid(True);
ax1.set_xlabel("x/c [-]");
ax1.set_ylabel("Cp [-]");

ax2 = plt.subplot(1, 2, 2);
ax2.plot(X, M_local, 'r');
ax2.set_xlabel("x/c [-]");
ax2.set_ylabel("M [-]");
ax2.grid(True);
plt.show();
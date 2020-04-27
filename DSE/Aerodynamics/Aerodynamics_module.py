from numpy import*
from matplotlib import pyplot as plt


def XFLR5PolarDataextraction(path_polars, path_X):
	polars = array(open(f_Polars).read().splitlines());
	unwantedtxt = [polars[0], polars[1], polars[2], polars[3], polars[4], polars[6], polars[7]];
	X = genfromtxt(f_X)[:, 0];

	for i in range(len(unwantedtxt)):
		idx = where(polars == unwantedtxt[i]);
		if i == len(unwantedtxt) - 1:
			Cp_start_idx = idx[0] - range(0, len(idx[0]));
		polars = delete(polars, idx, 0);

	polars = ([x.split() for x in polars]);
	polars = array([[float(j) for j in i] for i in (polars)]);
	return polars, X, Cp_start_idx;

def ApplyCompressibilityCorrection(M, C):
	beta = (sqrt(1 - M**2));
	coeff = 1/(beta)	# + (M**2/(1 + beta))*C/2);
	print(amax(coeff));
	return C*coeff;

def ComputeLocalPressure(M, Cp, p_inf, gamma, ratio):
	P = ((Cp*gamma*M**2)/2 + 1)*p_inf;
	P_tot = p_inf*(1 + (gamma - 1)/2*M**2)**(gamma/(gamma - 1));
	print(P_tot);
	if ratio == True:
		return P_tot/P;
	elif ratio == False:
		return P;

def IsentropicRelation(gamma, P_or_T_ratio, ratio_value):
	if P_or_T_ratio == "Pressure":
		ratio_value = sign(ratio_value)*absolute(ratio_value)**((gamma - 1)/gamma);
	M = sqrt((ratio_value - 1)*2/(gamma - 1));
	return M;

#%% Input data
f_Polars = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\Polars2.txt";
f_X = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\x_col.txt";
M = 0.78;
p_inf = 100000.;
gamma = 1.4;

polars, X, Cp_start_idx = XFLR5PolarDataextraction(f_Polars, f_X);
data = vstack((polars)[Cp_start_idx - 1]);
AOA = data[:, 0];
Cl = data[:, 2];	Cl_comp = ApplyCompressibilityCorrection(M, Cl);
Cd = data[:, 1];	Cd_comp = ApplyCompressibilityCorrection(M, Cd);
Cm = data[:, 3];	Cm_comp = ApplyCompressibilityCorrection(M, Cm);

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
ax1.set_xlabel("alpha [degrees]")
ax1.set_ylabel("Cl [-]")

ax2 = plt.subplot(2, 2, 2);
ax2.plot(AOA, Cd, label = "Incompressible");
ax2.plot(AOA, Cd_comp, 'r', label = "Compressible M =" + str(M));
ax2.set_xlabel("alpha [degrees]")
ax2.set_ylabel("Cd [-]")

ax3 = plt.subplot(2, 2, 3);
ax3.plot(AOA, Cm, label = "Incompressible");
ax3.plot(AOA, Cm_comp, 'r', label = "Compressible M =" + str(M));
ax3.set_xlabel("alpha [degrees]")
ax3.set_ylabel("Cm [-]")

ax4 = plt.subplot(2, 2, 4);
ax4.plot(Cd, Cl, label = "Incompressible");
ax4.plot(Cd_comp, Cl_comp, 'r', label = "Compressible M =" + str(M));
ax4.set_xlabel("Cd [-]")
ax4.set_ylabel("Cl [-]")
ax1.grid(True);
ax2.grid(True);
ax3.grid(True);
ax4.grid(True);
ax1.legend();
ax2.legend();
ax3.legend();
ax4.legend();

aoa = int(input("Specify Angle of attack :"));
req_idx = where(AOA == aoa);
cp = CP[:, req_idx].reshape(1, -1)[0];
cp_comp = ApplyCompressibilityCorrection(M, cp);
P_local = ComputeLocalPressure(M, cp_comp, p_inf, gamma, True);
M_local = IsentropicRelation(gamma, "Pressure", P_local);

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
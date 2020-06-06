from numpy import*
from matplotlib import pyplot as plt
from numpy.linalg import inv, det
from Aerodynamics.Input_parameters import*

#%% ------------------- Functions ----------------------
def Lagrange_Basis(nodes, f, plot_nodes):
	""" This function interpolates the function values "f"
	at the given nodes "nodes", and evaluates the interpolant
	at the plot_nodes
	Input:
		nodes = interpolation nodes
		f = function value at interpolation nodes
		plot_nodes = nodes at which the interpolant is evaluated
	Output:
		F = The interpolant evaluated at "plot_nodes"
	"""
	N = len(nodes);
	#basis = zeros((N), dtype = str);
	basis = [];
	F = 0;
	for i in range(N):
		l_i = 1;
		for j in range(N):
			if i != j:
				l_i *= (plot_nodes - nodes[j])/(nodes[i] - nodes[j]);
				basis.append("(x - {})/{} * ".format(nodes[j], nodes[i] - nodes[j]));
		F += l_i*f[i];
		if i != N - 1: basis.append(str(f[i]) + " + ");
		else: basis.append(str(f[i]));
	#for i in range(len(basis)):
	#	print(basis[i], end = '');
	#	if i % 3 == 0: print("\n");
	return F;

def RBF(r):
	""" Defining the Radial Basis Function
	Input:
		r = nodal distance
	Output:
		phi = RBF evaluated at r
	"""
	epi = 260;
	phi = sqrt(1 + (epi*r)**2);
	return phi;

def RBF_1DInterpol(nodes, f, plot_nodes):
	""" Function to Generate Matrix for RBF interpolation
	Input Arguments:
		nodes = interpolation nodes (numpy array)
	Output:
		F = The interpolant evaluated at "plot_nodes"
	"""
	X, XT = meshgrid(nodes, nodes);
	r = XT - X;
	A = RBF(r);
	coeff = inv(A).dot(f);

	X_p, n_p = meshgrid(plot_nodes, nodes);
	F_mat = RBF(X_p - n_p)*coeff.reshape(len(coeff), 1);
	F = sum(F_mat, axis = 0);
	return F;




#%% --------------- Main -----------------------
dir_CL = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Structures\liftdistribution.txt";
dir_Cm = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Structures\troquedistribution.txt";
data_CL = genfromtxt(dir_CL);
data_Cm = genfromtxt(dir_Cm);
CL = data_CL[1, :];												# CL values
Cm = data_Cm[1, :];												# Cm values around centre of gravity
span_location = data_CL[0, :];									# Spanwise stations
Lift_dist = CL*0.5*rho_cruise*V_cruise**2*ave_chord;			# Distributed lift load in N/m

x = linspace(span_location[0], span_location[-1], 100);
CL_lag = Lagrange_Basis(span_location, Lift_dist, x);
CL_rbf = RBF_1DInterpol(span_location, Lift_dist, x);
plt.figure(figsize = (18, 8));
plt.plot(span_location, Lift_dist, "x");
plt.plot(x, CL_lag, label = "Lagrange Interpolation");
plt.plot(x, CL_rbf, label = "RBF Interpolation");
plt.grid(True);
plt.legend();
plt.show();


from numpy import*
from matplotlib import pyplot as plt
from numpy.linalg import inv
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
	F = 0;
	for i in range(N):
		l_i = 1;
		for j in range(N):
			if i != j:
				l_i *= (plot_nodes - nodes[j])/(nodes[i] - nodes[j]);
		F += l_i*f[i];
	coeff = f;
	return F, coeff;

def RBF(r):
	""" Defining the Radial Basis Function
	Input:
		r = nodal distance
	Output:
		phi = RBF evaluated at r
	"""
	epi = 10;
	phi = sqrt(1 + (epi*r)**2);
	return phi;

def RBF_1DInterpol(nodes, f, plot_nodes, reconstruct = None):
	""" Function to Generate Matrix for RBF interpolation
	Input Arguments:
		nodes = interpolation nodes (numpy array)
	Output:
		F = The interpolant evaluated at "plot_nodes"
	"""
	if isinstance(reconstruct, ndarray):
		coeff = reconstruct;
	else:
		X, XT = meshgrid(nodes, nodes);
		r = XT - X;
		A = RBF(r);
		coeff = inv(A).dot(f);

	X_p, n_p = meshgrid(plot_nodes, nodes);
	F_mat = RBF(X_p - n_p)*coeff.reshape(len(coeff), 1);
	F = sum(F_mat, axis = 0);
	return F, coeff;

## Quadrature
def Quadrature_weights(x):
	""" Function to compute the weights for the Quadrature rule
	Input Arguments:
		x = Quadrature nodes (numpy float array)
	Output:
		w = Weights of the Quadrature rule
	"""
	N = len(x);
	A = ones((N, N));
	B = zeros((N, 1));
	for i in range(1, N):
		A[i, :] = x**i;
	
	a = x[0];		b = x[-1];
	for i in range(1, N + 1):
		B[i - 1] = (b**i - a**i)/i;
	w = inv(A).dot(B);
	return x, w.T[0];

def Transform_Quadrature(x, w, a, b):
	x_t = (b - a)/2*(x + 1) + a;
	w_t = (b - a)/2*w;
	return x_t, w_t;

def Quadrature_Integral(f, w, dim):
	""" Function to Compute 1D and 2D integral using Quadrature rule
	Input Arguments:
		f = Matrix containing function values at the nodes (numpy array)
		w = Quadrature weights, if dim = 2D then w has to have 2 rows (numpy array)
		dim = Dimention of definite integral, 1D or 2D (str)
	Output:
		I = Value of Integral
	"""
	if dim == "1D":
		I = f.dot(w);
	elif dim == "2D":
		I = f.dot(w[0, :]).dot(w[1, :]);
	return I;

def Generate_MomentShear_Diagram(F, F_args, y0, y1, DOP):
	y_intg_std = linspace(-1, 1, DOP + 1);
	_, w_intg_std = Quadrature_weights(y_intg_std);

	y = linspace(y0, y1, DOP + 1);
	V = zeros_like(y);
	M = zeros_like(y);
	for i in range(len(y)):
		y_out, w_out = Transform_Quadrature(y_intg_std, w_intg_std, y[len(y) - 1 - i], y1);
		g = zeros_like(y);
		for j in range(len(y_out)):
			eta_in, w_in = Transform_Quadrature(y_intg_std, w_intg_std, y_out[len(y_out) - 1 - j], y1);
			#F_args = (F_args[0], F_args[1], eta_in, F_args[2]);
			g[j] = Quadrature_Integral(-F(F_args[0], F_args[1], eta_in, F_args[2])[0], w_in, "1D");
		V[i] = g[-1];
		M[i] = Quadrature_Integral(-g, w_out, "1D");
	return M, V, y[::-1];


#%% --------------- Main -----------------------
dir_CL = "liftdistribution.txt";
dir_Cm = "troquedistribution.txt";
data_CL = genfromtxt(dir_CL);
data_Cm = genfromtxt(dir_Cm);
CL = data_CL[1, :];													# CL values
Cm = data_Cm[1, :];													# Cm values around centre of gravity
centre_pressure = data_Cm[2, :];									# x location of centre of pressure measured from LE root chord
span_location = data_CL[0, :];										# Spanwise stations
#FF = lambda x: 1 - 1/span_location[-1]*x;							# Verification
#Lift_dist = FF(span_location);										# Distributed lift load in N/m for verification
Lift_dist = CL*0.5*rho_cruise*V_cruise**2;							# Distributed lift load in N/m
Moment_dist = Cm*0.5*rho_cruise*V_cruise**2*ave_chord;				# Distributed pitching moment load in Nm/m

y = linspace(span_location[0], span_location[-1], 100);				# plot_nodes
CL_lag, _ = Lagrange_Basis(span_location, Lift_dist, y);			# Interpolate Distributed lift with Lagrange basis functions
CL_rbf, coeff = RBF_1DInterpol(span_location, Lift_dist, y);		# Interpolate Distributed lift with radial basis functions
Cm_lag, _ = Lagrange_Basis(span_location, Moment_dist, y);			# Interpolate Distributed pitching moment with Lagrange basis functions
Cm_rbf, _ = RBF_1DInterpol(span_location, Moment_dist, y);			# Interpolate Distributed pitching moment with radial basis functions

y_half = span_location[int(len(span_location)/2):len(span_location)];
M, V, y_m = Generate_MomentShear_Diagram(RBF_1DInterpol, [span_location, Lift_dist, coeff], y_half[0], y_half[-1], 11);

plt.figure(figsize = (18, 8));
plt.plot(span_location, Lift_dist, "x");
plt.plot(y, CL_lag, label = "Lagrange Interpolation");
plt.plot(y, CL_rbf, label = "RBF Interpolation");
plt.ylabel("Distributed Lift Load [N/m]");
plt.xlabel("Span [m]");
plt.grid(True);
plt.legend();

#plt.figure(figsize = (18, 8));
#plt.plot(span_location, Moment_dist, "x");
#plt.plot(y, Cm_lag, label = "Lagrange Interpolation");
#plt.plot(y, Cm_rbf, label = "RBF Interpolation");
#plt.ylabel("Distributed Pitching Moment Load [N/m]");
#plt.xlabel("Span [m]");
#plt.grid(True);
#plt.legend();

plt.figure(figsize = (18, 8));
plt.plot(y_m, M, label = "Internal Moment [Nm]");
plt.plot(y_m, V, label = "Shear [N]");
plt.ylabel("Interal Load Distributed");
plt.xlabel("y [m]");
plt.grid(True);
plt.legend();


plt.show();


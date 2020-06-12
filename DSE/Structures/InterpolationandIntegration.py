from numpy import*
from matplotlib import pyplot as plt
from numpy.linalg import inv
from Aerodynamics.Input_parm import*

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

def Generate_MomentShear_Diagram_Dist_Lift(F, F_args, y0, y1, DOP):
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

def Generate_MomentShear_Diagram_Engine_strut(M_lift, V_lift, y, T_e, W_e, LE_sweep, degree_strut, y_e, z_e, y_s, z_s):
	L = -min(V_lift);
	M_L = max(M_lift);

	x_e = y_e*tan(LE_sweep) - 1.768;									# engine distance from root start
	x_s = y_s*tan(LE_sweep);											#

	#2D FBD based on FBD in discord channel
	#S_z = (W_e*y_e - M_L)/((z_s/tan(degree_strut)) - y_s);
	#S_y = S_z/tan(degree_strut);
	S_z = -(W_e*y_e - M_L)/y_s;
	S_y = S_z/tan(degree_strut);
	A_z = L - W_e - S_z;
	A_y = -S_y;
	print(S_y, S_z);
	M_e = S_z*(y_s - y_e);												# moment at engine
	M_engine = W_e*y_e;
	M_r = S_z*(y_s) + M_engine;											# moment at root, should be equal to M_L, but isn't

	V_new = zeros(len(y_half));
	M_new = zeros(len(y_half));

	for i in range(len(y_half)):
		if y_half[i] == 0:
			V_new[i] = S_z + W_e;										# wing loading check
			M_new[i] = -M_r;											# moment check, doesn't work
		elif y_half[i] < y_e and y_half[i] > 0:
			V_new[i] = S_z + W_e;
			M_new[i] = (M_engine)/y_e*(y_half[i] - y_e) + (M_e)/(y_s - y_e)*(y_half[i] - y_s);
		elif y_half[i] < y_s and y_half[i] > y_e:
			V_new[i] = S_z;
			M_new[i] = (M_e)/(y_s - y_e)*(y_half[i] - y_s);

	V = V_new + V_lift;													# adding wing, strut and engine shear
	M = M_new + M_lift;													# adding wing, strut and engine moments #does not work
	return V_new, M_new, V, M;

def Generate_Torque_diagram(torque_dist_f, f_args, y0, y1, DOP):
	y_intg_std = linspace(-1, 1, DOP + 1);
	_, w_intg_std = Quadrature_weights(y_intg_std);

	y = linspace(y0, y1, DOP + 1);
	T = zeros_like(y);
	for i in range(len(y)):
		y_out, w_out = Transform_Quadrature(y_intg_std, w_intg_std, y[len(y) - 1 - i], y1);
		g = zeros_like(y);
		for j in range(len(y_out)):
			eta_in, w_in = Transform_Quadrature(y_intg_std, w_intg_std, y_out[len(y_out) - 1 - j], y1);
			g[j] = Quadrature_Integral(-torque_dist_f(f_args[0], f_args[1], eta_in, f_args[2])[0], w_in, "1D");
		T[i] = g[-1];
	return T, y[::-1];


#%% ----------------- Main -----------------------
dir_CL = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Structures\liftdistribution.txt";			# Chnage to your path
dir_Cm = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Structures\troquedistribution.txt";
data_CL = genfromtxt(dir_CL);
data_Cm = genfromtxt(dir_Cm);
CL = data_CL[1, :];													# CL values
Cm = data_Cm[1, :];													# Cm values around centre of gravity
centre_pressure = data_Cm[2, :];									# x location of centre of pressure measured from LE root chord
span_location = data_CL[0, :];										# Spanwise stations
#FF = lambda x: 1 - 1/span_location[-1]*x;							# Verification
#Lift_dist = FF(span_location)*1000;								# Distributed lift load in N/m for verification
Lift_dist = CL*0.5*rho_cruise*V_cruise**2;							# Distributed lift load in N/m
Moment_dist = Cm*0.5*rho_cruise*V_cruise**2*ave_chord;				# Distributed pitching moment load in Nm/m


y = linspace(span_location[0], span_location[-1], 100);				# plot_nodes
CL_lag, _ = Lagrange_Basis(span_location, Lift_dist, y);			# Interpolate Distributed lift with Lagrange basis functions
CL_rbf, coeff = RBF_1DInterpol(span_location, Lift_dist, y);		# Interpolate Distributed lift with radial basis functions
Cm_lag, _ = Lagrange_Basis(span_location, Moment_dist, y);			# Interpolate Distributed pitching moment with Lagrange basis functions
Cm_rbf, _ = RBF_1DInterpol(span_location, Moment_dist, y);			# Interpolate Distributed pitching moment with radial basis functions
centre_pressure, coeff_cp = RBF_1DInterpol(span_location, centre_pressure, y);

y_half = linspace(0, span_location[-1], 100);
M, V, y_m = Generate_MomentShear_Diagram_Dist_Lift(RBF_1DInterpol, [span_location, Lift_dist, coeff], y_half[0], y_half[-1], 9);

#%% ------------------- FBD Solve -------------------
V_lift, _ = RBF_1DInterpol(y_m, V, y_half);
M_lift, _ = RBF_1DInterpol(y_m, M, y_half);
T_e = 106000.;
W_e = 2000*9.81;
LE_sweep = radians(30);
y_e = 4.893;														# engine distance from fuselage
z_e = 1.434;														# engine thrust location under wing
y_s = 17.720;														# strut location
z_s = 0.732;														# length of vertical part of strut
degree_strut = radians(13.26);										# strut vertical angle
V_new, M_new, V, M = Generate_MomentShear_Diagram_Engine_strut(M_lift, V_lift, y_half, T_e, W_e, LE_sweep, degree_strut, y_e, z_e, y_s, z_s);


#%% ------------------ Plotting ------------------------
#plt.figure(figsize = (18, 8));
#plt.plot(span_location, Lift_dist, "x");
#plt.plot(y, CL_lag, label = "Lagrange Interpolation");
#plt.plot(y, CL_rbf, label = "RBF Interpolation");
#plt.ylabel("Distributed Lift Load [N/m]");
#plt.xlabel("Span [m]");
#plt.grid(True);
#plt.legend();

#plt.figure(figsize = (8, 5));
#plt.plot(y_half, M_lift, label = "Internal Moment due to lift [Nm]");
#plt.plot(y_half, V_lift, label = "Internal Shear due to lift [N]");
#plt.ylabel("Interal Load Distributed");
#plt.xlabel("y [m]");
#plt.grid(True);
#plt.legend();


#plt.figure(figsize = (8, 5));
#plt.plot(y_half, M_new, label = "Internal Moment due to engine and strut [Nm]");
#plt.plot(y_half, V_new, label = "Internal Shear due to engine and strut [N]");
#plt.ylabel("Interal Load Distributed");
#plt.xlabel("y [m]");
#plt.grid(True);
#plt.legend();

#plt.figure(figsize = (18, 8));
#plt.plot(y_half, M, label = "Internal Moment (total) [Nm]");
#plt.plot(y_half, V, label = "Internal Shear (total) [N]");
#plt.ylabel("Interal Load Distributed");
#plt.xlabel("y [m]");
#plt.grid(True);
#plt.legend();

#plt.show();

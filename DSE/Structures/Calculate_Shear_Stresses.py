from numpy import*
from InterpolationandIntegration import*
from Wingbox_MOI import*
from Aerodynamics.Aero_tools import sweep_x
from matplotlib import pyplot as plt

def Compute_sect_maxShear(V, T, tau_allow, t_c, c_arr, t_arr, n_stif_top, n_stif_bot, A_stif, A_spar_cap, iter = False):
	if iter == True:
		print("Itertating thickness");
		thickness = t_arr;
		t_min = zeros(len(c_arr));
		for chord in range(len(c_arr)):
			h = t_c*c_arr[chord];
			c = c_arr[chord];
			t_sec = [];
			tau = [];
			for t1 in thickness:
				for t2 in thickness:
					for t3 in thickness:
						centroid = calc_centroid(c, h, t1, t3, t2, 0, 0, 0, 0);
						#nodes = array([0, (h/2 - centroid), h]);
						#_, w = Quadrature_weights(nodes);
						Ixx, _ = calc_MOI(array([c]), array([h]), t1, t3, t2, array([0]), array([0]), 0, 0, centroid);
						q1 = -V[chord]/Ixx*(-t1*(0.5*h - centroid)*0.5*c);
						q2 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)**2 + ((0.5*h - centroid)**2)/2));
						q3 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)*h + (h**2)/2));
						#print((2*Quadrature_Integral(array([q1[0], q2[0], q3[0]]), w, "1D") - V[chord])/V[chord]);				# Quadrature integral for verfication
						q1 -= (T[chord]/(2*c*h));		q2 -= (T[chord]/(2*c*h));		q3 -= (T[chord]/(2*c*h));

						tau1 = abs(q1/t1);		tau2 = abs(q2/t2);		tau3 = abs(q3/t3);
						if tau1 < tau_allow and tau2 < tau_allow and tau3 < tau_allow:
							t_sec.append(max(t1, t2, t3));
							tau.append(min(tau1, tau2, tau3));
			t_min[chord] = min(t_sec);
			print("For section %d minimum required t = %f mm for stress of %f MPa" %(chord, t_min[chord]*1000, tau[chord]/1e6));
		return t_min;

	elif iter == False:
		t_top = t_arr[0, :];
		t_bot = t_arr[1, :];
		t_spar = t_arr[2, :];
		h_arr = t_c*c_arr;
		centroid = calc_centroid(c_arr, h_arr, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap);
		Ixx, _ = calc_MOI(c_arr, h_arr, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap, centroid);
		tau = zeros((len(c_arr), 4));
		for i in range(len(c_arr)):
			c = c_arr[i];		h = h_arr[i];
			t1 = t_top[i];		t2 = t_spar[i];			t3 = t_bot[i];
			q1 = -V[i]/Ixx[i]*(-t1*(0.5*h - centroid[i])*0.5*c);
			q2 = q1 - V[i]/Ixx[i]*(t2*(-(0.5*h - centroid[i])**2 + ((0.5*h - centroid[i])**2)/2));
			q3 = q1 - V[i]/Ixx[i]*(t2*(-(0.5*h - centroid[i])*h + (h**2)/2));
			q1 = abs(q1) + abs(T[i]/(2*c*h));		q2 = abs(q2) + abs(T[i]/(2*c*h));		q3 = abs(q3) + abs(T[i]/(2*c*h));
			tau1 = abs(q1/t1);		tau2 = abs(q2/t2);		tau3 = abs(q3/t3);
			tau[i, :] = asarray([max(tau1, tau2, tau3), tau1, tau3, tau2]);
		return tau;


def chord_y(c_root, c_tip, span, LE_sweep, y):
	f1 = tan(LE_sweep)*y;
	TE_sweep = sweep_x(1, LE_sweep, c_root, c_tip, span);
	f2 = c_root + tan(TE_sweep)*y;
	chord = f2 - f1;
	return chord;


#%% -------------- Main Inputs For iter == True ----------------
#c_root = 4.7;
#c_tip = 1.7;
#span = 50;
#LE_sweep = pi/6;
#LE_frac_unuse = 0.1;		TE_frac_unuse = 0.25;
#c_frac = 1 - LE_frac_unuse - TE_frac_unuse;
#no_stations = 10;
#y = linspace(0, y_half[-1], no_stations);
#V_, coeff_V = RBF_1DInterpol(y_half, V, y);
#chord = chord_y(c_root, c_tip, span, LE_sweep, y)*c_frac;
#CL, _ = RBF_1DInterpol(span_location, Lift_dist, y, reconstruct = coeff);
#centre_pressure = y*tan(LE_sweep) + 0.25*chord;
#SC = y*tan(LE_sweep) + (LE_frac_unuse + c_frac/2)*chord;
#torque_load = CL*(SC - centre_pressure);								# Distributed torque Nm/m
#torque_load_f, coeff_T = RBF_1DInterpol(y, torque_load, y);
#internal_torque, y_t = Generate_Torque_diagram(RBF_1DInterpol, [y, torque_load_f, coeff_T], y[0], y[-1], 9);
#internal_torque = internal_torque[::-1];
#t_arr = linspace(0.001, 0.02, 20);
#t_min_lst = Compute_sect_maxShear(V_, internal_torque, 200e6, 0.14*0.8, chord, t_arr, 0, 0, 0, 0, iter = True);

##%%
#fig = plt.figure(figsize = (15, 5));
#ax1 = plt.subplot(1, 3, 1);
#ax1.plot(y, t_min_lst*1000);
#ax1.grid(True);
#ax1.set_xlabel("Span location [m]");
#ax1.set_ylabel("Spar thickness [mm]");

#ax2 = plt.subplot(1, 3, 2);
#ax2.plot(y, V_/1000);
#ax2.grid(True);
#ax2.set_xlabel("Span location [m]");
#ax2.set_ylabel("Shear force [kN]");

#ax3 = plt.subplot(1, 3, 3);
#ax3.plot(y, internal_torque/1000);
#ax3.grid(True);
#ax3.set_xlabel("Span location [m]");
#ax3.set_ylabel("Torque [kNm]");
#plt.show();

#%% -------------- Main Inputs For iter == False ----------------
#c_root = 4.7;																# Root chord
#c_tip = 1.7;																# Tip chord
#span = 50;																	# Span
#LE_sweep = pi/6;															# Leading edge sweep
#LE_frac_unuse = 0.1;		TE_frac_unuse = 0.25;							# Avalable chord fraction avaliable for the wing box
#c_frac = 1 - LE_frac_unuse - TE_frac_unuse;
#no_stations = 10;															# Number of spanwise stations
#y = linspace(0, y_half[-1], no_stations);									# Span locations
#V_, _ = RBF_1DInterpol(y_half, V, y);										# Shear force at the span locations
#chord = chord_y(c_root, c_tip, span, LE_sweep, y)*c_frac;					# Chord length at the span locations
#CL, _ = RBF_1DInterpol(span_location, Lift_dist, y, reconstruct = coeff);	# Total lift distribution
#centre_pressure = y*tan(LE_sweep) + 0.25*chord;								# Centre of pressure
#SC = y*tan(LE_sweep) + (LE_frac_unuse + c_frac/2)*chord;					# Location of shear centre, assuming midpoint of box
#torque_load = CL*(SC - centre_pressure);									# Distributed torque Nm/m
#torque_load_f, coeff_T = RBF_1DInterpol(y, torque_load, y);
#internal_torque, y_t = Generate_Torque_diagram(RBF_1DInterpol, [y, torque_load_f, coeff_T], y[0], y[-1], 9);
#internal_torque = internal_torque[::-1];									# Internal torque distribution

#t_top = np.array([0.0005,0.0005,0.001,0.0015,0.002,0.0015,0.0005,0.0025,0.0015,0.0010])
#t_bot = np.array([0.001,0.0025,0.0025,0.0025,0.0025,0.0025,0.0015,0.0015,0.0005,0.0005])
#t_spar = np.array([0.002,0.003,0.004,0.005,0.005,0.004,0.003,0.002,0.001,0.001])
#n_stif_top = np.array([4,4,4,4,4,4,4,21,17,11])
#n_stif_bot = np.array([4,16,19,20,19,17,20,10,4,2])
#A_stif = 0.002*0.05*0.04
#A_spar_cap = 0.002*0.07*0.05
#tau_max = Compute_sect_maxShear(V_, internal_torque, 200e6, 0.14*0.8, chord, vstack((t_top, t_bot, t_spar)), n_stif_top, n_stif_bot, A_stif, A_spar_cap, iter = False);


#%%
#fig = plt.figure(figsize = (15, 5));
#ax1 = plt.subplot(1, 3, 1);
#ax1.plot(y, tau_max/1e6);
#ax1.grid(True);
#ax1.set_xlabel("Span location [m]");
#ax1.set_ylabel("max Shear Stree [MPa]");

#ax2 = plt.subplot(1, 3, 2);
#ax2.plot(y, V_/1000);
#ax2.grid(True);
#ax2.set_xlabel("Span location [m]");
#ax2.set_ylabel("Shear force [kN]");

#ax3 = plt.subplot(1, 3, 3);
#ax3.plot(y, internal_torque/1000);
#ax3.grid(True);
#ax3.set_xlabel("Span location [m]");
#ax3.set_ylabel("Torque [kNm]");
#plt.show();






#%% Verification
#span = 2.5;
#y = linspace(0, span, 100);
#chord = ones_like(y);
#torque_load_f, coeff_T = RBF_1DInterpol(y, 20*1000*ones_like(y), y);
#internal_torque, y_t = Generate_Torque_diagram(RBF_1DInterpol, [y, torque_load_f, coeff_T], y[0], y[-1], 9);

#t_min_lst = Compute_sect_maxShear(zeros_like(y), internal_torque, 200e6, 0.25, chord);

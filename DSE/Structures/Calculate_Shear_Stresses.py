from numpy import*
from InterpolationandIntegration import*
from Wingbox_MOI import*
from Aerodynamics.Aero_tools import sweep_x
from matplotlib import pyplot as plt

def Compute_sect_maxShear(V, T, tau_allow, t_c, c_arr):
	thickness = linspace(0.0001, 0.02, 20);
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
					nodes = array([0, (h/2 - centroid), h]);
					_, w = Quadrature_weights(nodes);
					Ixx, _ = calc_MOI(array([c]), array([h]), t1, t3, t2, array([0]), array([0]), 0, 0, centroid);
					q1 = -V[chord]/Ixx*(-t1*(0.5*h - centroid)*0.5*c);
					q2 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)**2 + ((0.5*h - centroid)**2)/2));
					q3 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)*h + (h**2)/2));
					#print((2*Quadrature_Integral(array([q1[0], q2[0], q3[0]]), w, "1D") - V[chord])/V[chord]);
					q1 = abs(q1) + abs(T[chord]/(2*c*h));		q2 = abs(q2) + abs(T[chord]/(2*c*h));		q3 = abs(q3) + abs(T[chord]/(2*c*h));
					
					tau1 = abs(q1/t1);		tau2 = abs(q2/t2);		tau3 = abs(q3/t3);
					if tau1 < tau_allow and tau2 < tau_allow and tau3 < tau_allow:
						t_sec.append(max(t1, t2, t3));
						tau.append(min(tau1, tau2, tau3));
		#print(t_sec)
		t_min[chord] = min(t_sec);
		print("For section %d minimum required t = %f mm for stress of %f MPa" %(chord, t_min[chord]*1000, tau[chord]/1e6));
	return t_min;

def chord_y(c_root, c_tip, span, LE_sweep, y):
	f1 = tan(LE_sweep)*y;
	TE_sweep = sweep_x(1, LE_sweep, c_root, c_tip, span);
	f2 = c_root + tan(TE_sweep)*y;
	chord = f2 - f1;
	return chord;




#%% -------------- Main ----------------
c_root = 4.7;
c_tip = 1.7;
span = 50;
LE_sweep = pi/6;
y = linspace(0, y_half[-1], 10);
V_, _ = RBF_1DInterpol(y_half, V, y);
chord = chord_y(c_root, c_tip, span, LE_sweep, y)*0.65;

CL, _ = RBF_1DInterpol(span_location, Lift_dist, y, reconstruct = coeff);
#centre_pressure, _ = RBF_1DInterpol(span_location, centre_pressure, y, reconstruct = coeff_cp);
centre_pressure = y*tan(LE_sweep) + 0.25*chord;
SC = y*tan(LE_sweep) + 0.425*chord;
torque_load = CL*(SC - centre_pressure);								# Distributed torque Nm/m
torque_load_f, coeff_T = RBF_1DInterpol(y, torque_load, y);
internal_torque, _ = Generate_Torque_diagram(RBF_1DInterpol, [y, torque_load_f, coeff_T], y[0], y[-1], 9);
internal_torque = internal_torque[::-1];

t_min_lst = Compute_sect_maxShear(V_, internal_torque, 200e6, 0.14*0.8, chord);
plt.plot(y, torque_load, "x");
plt.plot(y, torque_load_f);


#%% Verification
#span = 2.5;
#y = linspace(0, span, 100);
#chord = ones_like(y);
#torque_load_f, coeff_T = RBF_1DInterpol(y, 20*1000*ones_like(y), y);
#internal_torque, y_t = Generate_Torque_diagram(RBF_1DInterpol, [y, torque_load_f, coeff_T], y[0], y[-1], 9);

#t_min_lst = Compute_sect_maxShear(zeros_like(y), internal_torque, 200e6, 0.25, chord);


#%%
fig = plt.figure(figsize = (15, 5));
ax1 = plt.subplot(1, 3, 1);
ax1.plot(y, t_min_lst*1000);
ax1.grid(True);
ax1.set_xlabel("Span location [m]");
ax1.set_ylabel("Spar thickness [mm]");

ax2 = plt.subplot(1, 3, 2);
ax2.plot(y, V_/1000);
ax2.grid(True);
ax2.set_xlabel("Span location [m]");
ax2.set_ylabel("Shear force [kN]");

ax3 = plt.subplot(1, 3, 3);
ax3.plot(y, internal_torque/1000);
ax3.grid(True);
ax3.set_xlabel("Span location [m]");
ax3.set_ylabel("Torque [kNm]");
plt.show();

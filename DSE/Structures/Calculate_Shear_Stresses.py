from numpy import*
from InterpolationandIntegration import*
from Wingbox_MOI import*
from Aerodynamics.Aero_tools import sweep_x
from matplotlib import pyplot as plt

def Compute_sect_maxShear(V, tau_allow, t_c, c_arr):
	thickness = linspace(0.0001, 0.002, 20);
	t_min = zeros(len(c_arr));
	for chord in range(len(c_arr)):
		h = 0.8*t_c*c_arr[chord];
		c = 0.65*c_arr[chord];
		t_sec = [];
		tau = [];
		nodes = array([0, h/2, h]);
		_, w = Quadrature_weights(nodes);
		for t1 in thickness:
			for t2 in thickness:
				for t3 in thickness:
					centroid = calc_centroid(c, h, t1, t3, t2, 0, 0, 0, 0);
					Ixx, _ = calc_MOI(array([c]), array([h]), t1, t3, t2, array([0]), array([0]), 0, 0, centroid);
					q1 = -V[chord]/Ixx*(-t1*(0.5*h - centroid)*0.5*c);
					q2 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)**2 + ((0.5*h - centroid)**2)/2));
					q3 = q1 - V[chord]/Ixx*(t2*(-(0.5*h - centroid)*h + (h**2)/2));
					#print((2*Quadrature_Integral(array([q1[0], q2[0], q3[0]]), w, "1D") - V[chord])/V[chord]);
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
chord = chord_y(c_root, c_tip, span, LE_sweep, y);
t_min_lst = Compute_sect_maxShear(V_, 200e6, 0.14, chord);


fig = plt.figure();
ax1 = plt.subplot(1, 2, 1);
ax1.plot(y, t_min_lst*1000);
ax1.grid(True);
ax1.set_xlabel("Span location [m]");
ax1.set_ylabel("Spar thickness [mm]");

ax2 = plt.subplot(1, 2, 2);
ax2.plot(y, V_/1000);
ax2.grid(True);
ax2.set_xlabel("Span location [m]");
ax2.set_ylabel("Shear force [kN]");
plt.show();

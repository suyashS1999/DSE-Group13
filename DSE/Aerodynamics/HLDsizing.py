from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Aero_tools import sweep_x

class HighLiftDevices():
	def __init__(self):
		self.TE_HLDs = dict();
		self.LE_HLDs = dict();

		self.TE_HLDs["PlainFlap"] = [20, 60, 0.9];
		self.TE_HLDs["SingleSlottedFlap"] = [20, 40, 1.3];
		self.TE_HLDs["SingleFowlerFlap"] = [15, 40, (lambda delta, c_f_c: 1.3*(1 + (6e-3*delta + 0.41)*c_f_c))];
		self.TE_HLDs["DoubleFowlerFlap"] = [10, 50, (lambda delta, c_f_c: 1.6*(1 + (9.6e-3*delta + 0.406)*c_f_c))];
		self.TE_HLDs["TripleFowlerFlap"] = [20, 40, (lambda delta, c_f_c: 1.9*(1 + (9.6e-3*delta + 0.406)*c_f_c))];

		self.LE_HLDs["FixedSlot"] = [0, 0, 0.2];
		self.LE_HLDs["LEFlap"] = [0, 0, 0.3];
		self.LE_HLDs["KrugerFlap"] = [0, 0, 0.3];
		self.LE_HLDs["Slat"] = [0, 0, 0.4];

	def Delta_CL_alpha0(self, delta_Clmax, Swf_S, lambda_hingeLn, delta_Clmax_args):
		if isinstance(delta_Clmax, float):
			Delta_CL = 0.9*delta_Clmax*Swf_S*cos(lambda_hingeLn);
		else:
			Delta_CL = 0.9*delta_Clmax(*delta_Clmax_args)*Swf_S*cos(lambda_hingeLn);
		delta_alpha0 = -15;
		Delta_alpha0 = delta_alpha0*Swf_S*cos(lambda_hingeLn);
		return Delta_CL, Delta_alpha0;


	def plot_Delta_CL_alpha0(self, LE_or_TE, c_f_c):
		if LE_or_TE == "TE":
			HLDs = self.TE_HLDs;
		elif LE_or_TE == "LE":
			HLDs = self.LE_HLDs;
		Swf_S = linspace(0.1, 1.0, 10);
		lambda_hingeLn = linspace(0, pi/2, 10);
		Swf_S_mesh, lambda_hingeLn_mesh = meshgrid(Swf_S, lambda_hingeLn);
		for hld in list(HLDs):
			fig = plt.figure(figsize = (16, 6));
			ax1 = fig.add_subplot(1, 2, 1, projection = "3d");
			ax2 = fig.add_subplot(1, 2, 2, projection = "3d");
			HLD = list(HLDs[hld]);
			delta_TO = HLD[0];
			delta_Land = HLD[1];
			delta_Clmax = HLD[2];
			Delta_CL, Delta_alpha0 = self.Delta_CL_alpha0(delta_Clmax, Swf_S_mesh, lambda_hingeLn_mesh, (delta_Land, c_f_c));
			
			ax1.plot_surface(Swf_S_mesh, lambda_hingeLn_mesh, Delta_CL, rstride = 1, cstride = 1, cmap = "viridis", edgecolor = "none");
			ax2.plot_surface(Swf_S_mesh, lambda_hingeLn_mesh, Delta_alpha0, rstride = 1, cstride = 1, cmap = "viridis", edgecolor = "none");
			ax1.set_xlabel("Swf/S [-]");
			ax1.set_ylabel("Hinge line sweep [rad]");
			ax1.set_zlabel("Delta CL [-]");
			ax1.set_title(hld);
			ax2.set_xlabel("Swf/S [-]");
			ax2.set_ylabel("Hinge line sweep [rad]");
			ax2.set_zlabel("Delta alpha0 [-]");
			ax2.set_title(hld);
		return 0;

	def Adjust_CL_alpha(self, ori_CL_alpha, ori_CL0, ori_CL_max, LE_HLD, TE_HLD, c_f_c, Swf_S_TE, Swf_S_LE, lambda_hingeLn_LE, lambda_hingeLn_TE):
		fig = plt.figure(figsize = (9, 8));
		LE_hld = list(self.LE_HLDs[LE_HLD]);
		TE_hld = list(self.TE_HLDs[TE_HLD]);
		delta_Clmax_LE = LE_hld[2];
		delta_Clmax_TE = TE_hld[2];
		delta_Land = TE_hld[1];
		delta_TO = TE_hld[0];
		Delta_CL_LE, Delta_alpha0_LE = self.Delta_CL_alpha0(delta_Clmax_LE, Swf_S_LE, lambda_hingeLn_LE, (delta_Land, c_f_c));
		Delta_CL_TE, Delta_alpha0_TE = self.Delta_CL_alpha0(delta_Clmax_TE, Swf_S_TE, lambda_hingeLn_TE, (delta_Land, c_f_c));
		print(Delta_CL_LE);
		print(Delta_CL_TE);
		print(Delta_alpha0_LE);
		print(Delta_alpha0_TE);
		if TE_HLD == "SingleFowlerFlap" or TE_HLD == "DoubleFowlerFlap" or TE_HLD == "TripleFowlerFlap":
			if Swf_S_TE != 0: CL_alpha = (1 + c_f_c)*ori_CL_alpha;
			else: CL_alpha = ori_CL_alpha;
		else:
			CL_alpha = ori_CL_alpha;
		CL0 = ori_CL0;
		CL_max = ori_CL_max + Delta_CL_LE + Delta_CL_TE;
		print("# ---------------------------- #");
		print("New CL max = {}".format(CL_max));
		print("# ---------------------------- #");
		alpha = linspace(-5.0, 10, 20);
		CL_ori = ori_CL_alpha*alpha + ori_CL0;
		CL_new = CL_alpha*(alpha + Delta_alpha0_LE - Delta_alpha0_TE) + CL0;
		plt.plot(alpha, CL_new, label = "New Lift curve");
		plt.plot(alpha, CL_ori, label = "Original Lift curve");
		plt.grid(True);
		plt.xlabel(r"$\alpha$ [deg]");
		plt.ylabel(r"$C_L$ [-]");
		plt.legend(loc = "upper left");
		plt.xlim([-5, 24]);
		plt.ylim([-0.1, 1.1*CL_max]);
		plt.show();
		return 0;

def chord_y(c_root, c_tip, span, LE_sweep, y):
	f1 = tan(LE_sweep)*y;
	TE_sweep = sweep_x(1, LE_sweep, c_root, c_tip, span);
	f2 = c_root + tan(TE_sweep)*y;
	chord = f2 - f1;
	return chord;

def ComputeHLD_Dim(b0, Swf_S, S, sweep_LE, sweep_TE, root_c, tip_c):
	chord_b0 = chord_y(root_c, tip_c, span, sweep_LE, b0);
	bf = (-2*chord_b0 + sqrt((2*chord_b0)**2 - 4*(tan(sweep_TE) - tan(sweep_LE))*(-Swf_S*S)))/(2*(tan(sweep_TE) - tan(sweep_LE)));
	chord_b1 = chord_y(root_c, tip_c, span, sweep_LE, bf + b0);
	S_f = bf*(chord_b0 + chord_b1);
	return S_f, bf;

#%% ------------------ Main ------------------
span = 50.64;								# Wing span
root_c = 4.135518182;						# Root chord
tip_c = 1.819628;							# Tip chord
S = span*(root_c + tip_c)/2;				# Wing area
sweep_LE = pi/6;							# LE sweep
sweep_TE = sweep_x(1, sweep_LE, root_c, tip_c, span);
fuselage_d = 4.2;							# Fuselage diameter
fuselage_clear_TE = 1.5;					# Where the TE hld starts
fuselage_clear_LE = 0.5;					# Where the LE hld starts
HLD = HighLiftDevices();
#HLD.plot_Delta_CL_alpha0("LE", 0.25);
#HLD.plot_Delta_CL_alpha0("TE", 0.25);
#plt.show();
ori_CL_alpha = 0.08772433;					# Original Clean CL-alpha
ori_CL0 = 0.280878;							# Original Clean CL at zero aoa
ori_CL_max = 1.7536;						# Original Clean CL max
LE_HLD = "KrugerFlap";						# Choose leading edge HLD. Choose between ["FixedSlot", "LEFlap", "KrugerFlap"]
TE_HLD = "DoubleFowlerFlap";				# Choose trailing edge HLD. Choose between ["PlainFlap", "SingleSlottedFlap", "SingleFowlerFlap", "DoubleFowlerFlap", "TripleFowlerFlap"]
c_f_c = 0.25;								# fraction of the chord that is the HLD
Swf_S_TE = 0.5;								# Flapped (TE) area ratio
Swf_S_LE = 0.62;							# Flapped (LE) area ratio
lambda_hingeLn_LE = sweep_LE;				# Sweep of hinge line LE
lambda_hingeLn_TE = sweep_TE;				# Sweep of hinge line TE
HLD.Adjust_CL_alpha(ori_CL_alpha, ori_CL0, ori_CL_max, LE_HLD, TE_HLD, c_f_c, Swf_S_TE, Swf_S_LE, lambda_hingeLn_LE, lambda_hingeLn_TE);

b0_TE = fuselage_d/2 + fuselage_clear_TE;
b0_LE = fuselage_d/2 + fuselage_clear_LE;
HLD_LE_S, bf_LE = ComputeHLD_Dim(b0_LE, Swf_S_LE, S, sweep_LE, sweep_TE, root_c, tip_c);
HLD_TE_S, bf_TE = ComputeHLD_Dim(b0_TE, Swf_S_TE, S, sweep_LE, sweep_TE, root_c, tip_c);
print("# -------- LE HLD ----------#");
print("Half span = %f m \nFlapped area = %f m^2 \nHLD start span = %f m from fuselage centerline \nHLD end span = %f m from fuselage centerline" %(bf_LE, HLD_LE_S, b0_LE, b0_LE + bf_LE));
print("\n# -------- TE HLD ----------#");
print("Half span = %f m \nFlapped area = %f m^2 \nHLD start span = %f m from fuselage centerline \nHLD end span = %f m from fuselage centerline" %(bf_TE, HLD_TE_S, b0_TE, b0_TE + bf_TE));





		

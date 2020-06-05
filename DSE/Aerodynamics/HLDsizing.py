from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class HighLiftDevices():
	def __init__(self):
		self.TE_HLDs = dict();
		self.LE_HLDs = dict();

		self.TE_HLDs["PlainFlap"] = [20, 60, 0.9];
		self.TE_HLDs["SingleSlottedFlap"] = [20, 40, 1.3];
		self.TE_HLDs["SingleFowlerFlap"] = [15, 40, (lambda delta, c_f_c: 1.3*(1 + (6e-3*delta + 0.41)*c_f_c))];
		self.TE_HLDs["DoubleFowlerFlap"] = [20, 50, (lambda delta, c_f_c: 1.6*(1 + (9.6e-3*delta + 0.406)*c_f_c))];
		self.TE_HLDs["TripleFowlerFlap"] = [20, 40, (lambda delta, c_f_c: 1.9*(1 + (9.6e-3*delta + 0.406)*c_f_c))];

		self.LE_HLDs["FixedSlot"] = [0, 0, 0.2];
		self.LE_HLDs["LEFlap"] = [0, 0, 0.3];
		self.LE_HLDs["KrugerFlap"] = [0, 0, 0.3];

	def Delta_CL_alpha0(self, delta_Clmax, Swf_S, lambda_hingeLn, delta_Clmax_args):
		if isinstance(delta_Clmax, float):
			Delta_CL = 0.9*delta_Clmax*Swf_S*cos(lambda_hingeLn);
		else:
			Delta_CL = 0.9*delta_Clmax(*delta_Clmax_args)*Swf_S*cos(lambda_hingeLn);
		delta_alpha0 = -15;
		Delta_alpha0 = delta_alpha0*pi/180*Swf_S*cos(lambda_hingeLn);
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
		fig = plt.figure(figsize = (8, 8));
		LE_hld = list(self.LE_HLDs[LE_HLD]);
		TE_hld = list(self.TE_HLDs[TE_HLD]);
		delta_Clmax_LE = LE_hld[2];
		delta_Clmax_TE = TE_hld[2];
		delta_Land = TE_hld[1];
		Delta_CL_LE, Delta_alpha0_LE = self.Delta_CL_alpha0(delta_Clmax_LE, Swf_S_LE, lambda_hingeLn_LE, (delta_Land, c_f_c));
		Delta_CL_TE, Delta_alpha0_TE = self.Delta_CL_alpha0(delta_Clmax_TE, Swf_S_TE, lambda_hingeLn_TE, (delta_Land, c_f_c));
		if TE_HLD == "SingleFowlerFlap" or TE_HLD == "DoubleFowlerFlap" or TE_HLD == "TripleFowlerFlap":
			if Swf_S_TE != 0: CL_alpha = (1 + c_f_c)*ori_CL_alpha;
			else: CL_alpha = ori_CL_alpha;
		else:
			CL_alpha = ori_CL_alpha;
		CL0 = ori_CL0 - Delta_alpha0_LE - Delta_alpha0_TE;
		CL_max = ori_CL_max + Delta_CL_LE + Delta_CL_TE;
		print("# ---------------------------- #");
		print("New CL max = {}".format(CL_max));
		print("# ---------------------------- #");
		alpha = linspace(-5.0, 10, 20);
		CL_ori = ori_CL_alpha*alpha + ori_CL0;
		CL_new = CL_alpha*alpha + CL0;
		plt.plot(alpha, CL_new, label = "New Lift curve");
		plt.plot(alpha, CL_ori, label = "Original Lift curve");
		plt.grid(True);
		plt.xlabel("alpha [deg]");
		plt.ylabel("CL [-]");
		plt.legend();
		plt.show();
		return 0;

#%% ------------------ Main ------------------
HLD = HighLiftDevices();
#HLD.plot_Delta_CL_alpha0("LE", 0.25);
#HLD.plot_Delta_CL_alpha0("TE", 0.25);
#plt.show();
ori_CL_alpha = 0.074417;					# Original Clean CL-alpha
ori_CL0 = 0.300238;							# Original Clean CL at zero aoa
ori_CL_max = 1.113;							# Original Clean CL max
LE_HLD = "KrugerFlap";						# Choose leading edge HLD. Choose between ["FixedSlot", "LEFlap", "KrugerFlap"]
TE_HLD = "SingleFowlerFlap";				# Choose trailing edge HLD. Choose between ["PlainFlap", "SingleSlottedFlap", "SingleFowlerFlap", "DoubleFowlerFlap", "TripleFowlerFlap"]
c_f_c = 0.25;								# fraction of the chord that is the HLD
Swf_S_TE = 0.6;								# Flapped (TE) area ratio
Swf_S_LE = 0.8;								# Flapped (LE) area ratio
lambda_hingeLn_LE = 0.7;					# Sweep of hinge line LE
lambda_hingeLn_TE = 0.6;					# Sweep of hinge line TE
HLD.Adjust_CL_alpha(ori_CL_alpha, ori_CL0, ori_CL_max, LE_HLD, TE_HLD, c_f_c, Swf_S_TE, Swf_S_LE, lambda_hingeLn_LE, lambda_hingeLn_TE);


		

from numpy import*
from matplotlib import pyplot as plt

class HighLiftDevices():
	def __init__(self):
		self.TE_HLDs = dict();
		self.LE_HLDs = dict();

		self.TE_HLDs["PlainFlap"] = [20, 60, 0.9];
		self.TE_HLDs["SingleSlottedFlap"] = [20, 40, 1.3];
		self.TE_HLDs["SingleFowlerFlap"] = [15, 40, (lambda delta, c_f_c: 1.3*(1 + (6e-3*delta + 0.41)*c_f_c))];
		self.TE_HLDs["DoubleFowlerFlap"] = [20, 50, (lambda delta, c_f_c: 1.6*(1 + (9.6e-3*delta + 0.406)*c_f_c))];
		self.TE_HLDs["TripleFowlerFlap"] = [20, 40, (lambda delta, c_f_c: 1.9*(1 + (9.6e-3*delta + 0.406)*c_f_c))];

		#self.LE_HLDs["SlottedLEFlap"] = {}

	def Delta_CL_alpha0(self, LE_or_TE, ave_c):
		if LE_or_TE == "TE":
			HLDs = self.TE_HLDs;
			Swf_S = linspace(0.1, 1.0, 10);
			lambda_hingeLn = linspace(0, 0.5, 10);

			c_f_c = 0.3;
			fig = plt.figure(figsize = (8, 6));
			for hld in list(HLDs):
				HLD = list(HLDs[hld]);
				delta_TO = HLD[0];
				delta_Land = HLD[1];
				delta_Clmax = HLD[2];
				if isinstance(delta_Clmax, float):
					Delta_CL = 0.9*delta_Clmax*Swf_S*cos(lambda_hingeLn);
				else:
					Delta_CL = 0.9*delta_Clmax(delta_Land, c_f_c)*Swf_S*cos(lambda_hingeLn);
				plt.plot(Swf_S, Delta_CL, label = hld);
			plt.grid(True);
			plt.xlabel("Swf/S [-]");
			plt.ylabel("Delta CL [-]");
			plt.legend();
			plt.show();
		return 0;

HLD = HighLiftDevices();
HLD.Delta_CL_alpha0("TE", 2.89543, 40.63*pi/180);



		

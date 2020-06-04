from numpy import*
from numpy.linalg import inv, norm
import glob, os
from matplotlib import pyplot as plt
import time

class ExtractData_OpenVSP():
	""" Object to extract data from the output files made by VSPAero
	Input:
		dir = Directory to search for the OpenVSP data files
	"""
	def __init__(self, dir):
		self.file_types = [".polar", ".lod"];
		self.names = [];
		self.Polar_Dict = dict();							# Global dictionary for all the polars
		self.Load_Dict = dict();							# Global dictionary for all the load files
		os.chdir(dir);
		print("Found files:");
		for file_type in self.file_types:
			for file in glob.glob("*{}".format(file_type)):
				name = "%s%s" %(file_type[1:], file[: -len(file_type)]);
				self.names.append(name);
				print(file, "\n");
				if file.endswith(self.file_types[0]):		# Sorting out the Polars files
					polars = genfromtxt(file, dtype = "str");
					polars_data = asarray(polars[1::], dtype = float);
					parms = polars[0];
					Dict = dict();							# Local dictionary for a single polar
					for i in range(len(parms)):
						Dict[parms[i]] = polars_data[:, i];
					self.Polar_Dict[name] = Dict;			# Add local dictionary to global dictionary
				elif file.endswith(self.file_types[1]):		# Sorting out the Load files
					self.Load_Dict[name] = array(open(file).read().splitlines());
		self.Sort_LoadDistribution();

	# ------------------- Functions -------------------
	def plot_Polars(self):
		""" Function to plot the polars from the data
		Output:
			fig = plot of the polars
		"""
		fig = plt.figure(figsize = (12, 8));
		for name in (self.names):
			if name.startswith(self.file_types[0][1:]):
				alpha = self.Polar_Dict[name]["AoA"];
				CL = self.Polar_Dict[name]["CL"];
				CDtot = self.Polar_Dict[name]["CDtot"];
				L_D = self.Polar_Dict[name]["L/D"];
				Cm = self.Polar_Dict[name]["CMy"];

				label = name[len(self.file_types[0][1:]) :];
				ax1 = plt.subplot(2, 2, 1);									ax2 = plt.subplot(2, 2, 2);
				ax1.plot(alpha, CL, label = label);							ax2.plot(alpha, Cm, label = label);
				ax1.set_xlabel("alpha [degrees]");							ax2.set_xlabel("alpha [degrees]");
				ax1.set_ylabel("CL [-]");									ax2.set_ylabel("Cm [-]");
				ax1.grid(True);												ax2.grid(True);

				ax3 = plt.subplot(2, 2, 3);									ax4 = plt.subplot(2, 2, 4);
				ax3.plot(CDtot, CL, label = label);							ax4.plot(alpha, L_D, label = label);
				ax3.set_xlabel("CD [-]");									ax4.set_xlabel("alpha [degrees]");
				ax3.set_ylabel("CL [-]");									ax4.set_ylabel("CL/CD [-]");
				ax3.grid(True);												ax4.grid(True);

				ax1.legend();
				ax2.legend();
				ax3.legend();
				ax4.legend();
		return 0;

	def Sort_LoadDistribution(self):
		""" Function to extract the Lift and Drag distribution data 
		"""
		self.subDict = dict();
		for name in (self.names):
			if name.startswith(self.file_types[1][1:]):
				Data = self.Load_Dict[name];
				sections = where(Data == Data[0])[0];
				for i in range(len(sections)):
					if i == len(sections) - 1: data = Data[sections[i] + 3: len(Data) - 3];
					else: data = Data[sections[i] + 3: sections[i + 1] - 3];
					count = 1;
					secDict = dict();
					for j in range(len(data)):
						d = (data[j].split());
						if len(d) <= 3 and len(d) != 0: secDict[d[0]] = float(d[1]);		count += 1;
						elif len(d) > 0 and d[0] == "Comp": break;
						elif len(d) > 3 and j != count:
							for k in range(len(d)):
								try: 
									secDict[data[count].split()[k]].append(float(d[k]));
								except: 
									secDict[data[count].split()[k]] = ([float(d[k])]);
					self.subDict[secDict["AoA_"]] = secDict;


	def plot_LoadDistribution(self):
		""" Function to plot the Lift and Drag distribution over 
		the wing from the data
		Output:
			fig = plot of the distribution
		"""
		for name in (self.names):
			if name.startswith(self.file_types[1][1:]):
				fig = plt.figure(figsize = (12, 10));
				for i in range(len(self.subDict)):
					alpha = float(list(self.subDict)[i]);
					CL = asarray(self.subDict[alpha]["Cl"])*asarray(self.subDict[alpha]["Chord"])/asarray(self.subDict[alpha]["Cref_"]);
					CD = asarray(self.subDict[alpha]["Cd"])*asarray(self.subDict[alpha]["Chord"])/asarray(self.subDict[alpha]["Cref_"]);
					span = asarray(self.subDict[alpha]["Yavg"]);
					sort_idx = argsort(span);
					span = span[sort_idx];		CL = CL[sort_idx];		CD = CD[sort_idx];
					ax1 = plt.subplot(2, 1, 1);
					ax1.title.set_text("Lift Distribution {}".format(name[len(self.file_types[1][1:]) :]));
					ax1.plot(span, CL, "x-", label = "alpha = {}".format(alpha));
					ax1.set_xlabel("span [m]");
					ax1.set_ylabel("Cl [-]");
					ax1.grid(True);

					ax2 = plt.subplot(2, 1, 2);
					ax2.title.set_text("Drag Distribution {}".format(name[len(self.file_types[1][1:]) :]));
					ax2.plot(span, CD, "^-", label = "alpha = {}".format(alpha));
					ax2.set_ylabel("Cd [-]");
					ax2.set_xlabel("span [m]");
					ax2.grid(True);

					ax1.legend(loc = "upper right");
					ax2.legend(loc = "upper right");
		return 0;


	def Cm_CL_alpha_calc(self, precision = 25):
		""" Function to compute the Cm-alpha and CL-alpha curves.
		It uses linear regression to compute the curves. The precision is input does
		not need to be changed
		"""
		least_sq_Matrix = lambda x: vstack((ones(shape(x)), x));
		for name in self.names:
			if name.startswith(self.file_types[0][1:]):
				alpha = self.Polar_Dict[name]["AoA"];
				CL = self.Polar_Dict[name]["CL"];
				Cm = self.Polar_Dict[name]["CMy"];
				idx = asarray(range(0, len(alpha), int(precision/len(alpha))));
				alpha_i = alpha[idx];
				CL_i = CL[idx];
				Cm_i = Cm[idx];
				LSM = least_sq_Matrix(alpha_i);
				A = LSM.dot(LSM.T);

				a_CL = inv(A).dot(LSM.dot(CL_i));
				R2_CL = 1 - norm(LSM.T.dot(a_CL) - CL_i);
				a_Cm = inv(A).dot(LSM.dot(Cm_i));
				R2_Cm = 1 - norm(LSM.T.dot(a_Cm) - Cm_i);
				print("Test Case:", name[len(self.file_types[0][1:]) :])
				print("CL(alpha) = %f + %f alpha" %(a_CL[0], a_CL[1]), "\n R^2 = {}".format(R2_CL));
				print("Cm(alpha) = %f + %f alpha" %(a_Cm[0], a_Cm[1]), "\n R^2 = {}".format(R2_Cm));
				print("\n")
		return 0;

#%% ------------------- Input data -------------------
dir = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\OpenVSPSimData";			# Change this to you local directory
#dir = r"C:\Users\miksw\Desktop\DSE\DSE-Group13\DSE\Aerodynamics\OpenVSPSimData";
#%% ------------------- Main -------------------
vsp_data = ExtractData_OpenVSP(dir);
vsp_data.plot_Polars();
vsp_data.plot_LoadDistribution();
vsp_data.Cm_CL_alpha_calc();
plt.show();




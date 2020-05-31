from numpy import*
import glob, os
from matplotlib import pyplot as plt

class ExtractData_OpenVSP():
	""" Object to extract data from the output files made by VSPAero
	Input:
		dir = Directory to search for the OpenVSP data files
	"""
	def __init__(self, dir):
		self.file_types = [".polar"];			#, ".lod"
		self.files = [];
		self.names = [];
		os.chdir(dir);
		print("Found files:");
		for file_type in self.file_types:
			for file in glob.glob("*{}".format(file_type)):
				self.files.append(file);
				name = "%s%s" %(file_type[1:], file[: -len(file_type)]);
				self.names.append(name);
				print(self.files[-1], "\n");
				vars(self)[name] = genfromtxt(file);

	def plot_Polars(self):
		""" Function to plot plot the polars from the data
		Output:
			fig = plot of the polars
		"""
		fig = plt.figure(figsize = (16, 16));
		for name in (self.names):
			if name.startswith(self.file_types[0][1:]):
				polars = vars(self)[name][1::];
				alpha = polars[:, 2];
				CL = polars[:, 4];
				CDtot = polars[:, 7];
				L_D = polars[:, 9];
				Cm = polars[:, 17];
				label = name[len(self.file_types[0][1:]) :];
				ax1 = plt.subplot(2, 2, 1);									ax2 = plt.subplot(2, 2, 2);
				ax1.plot(alpha, CL, label = label);							ax2.plot(alpha, Cm, label = label);
				ax1.set_xlabel("alpha [degrees]");							ax2.set_xlabel("alpha [degrees]");
				ax1.set_ylabel("Cl [-]");									ax2.set_ylabel("Cm [-]");
				ax1.grid(True);												ax2.grid(True);

				ax3 = plt.subplot(2, 2, 3);									ax4 = plt.subplot(2, 2, 4);
				ax3.plot(CDtot, CL, label = label);							ax4.plot(alpha, L_D, label = label);
				ax3.set_xlabel("CD [-]");									ax4.set_xlabel("alpha [degrees]");
				ax3.set_ylabel("CL [-]");									ax4.set_ylabel("Cl/CD [-]");
				ax3.grid(True);												ax4.grid(True);

				ax1.legend();
				ax2.legend();
				ax3.legend();
				ax4.legend();
		plt.show();


#%% Input data
dir = r"C:\Users\Gebruiker\source\repos\DSE\DSE\Aerodynamics\OpenVSPSimData";			# Change this to you local directory
vsp_data = ExtractData_OpenVSP(dir);
vsp_data.plot_Polars();


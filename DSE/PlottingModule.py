from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm

def Plot_1D(X, Y, label = None, xlabel = "You forgot to label your axis", ylabel = "You forgot to label your axis", show = False, figure = None, xlim = None, ylim = None, figsize = None):
	""" This function will plot X, Y in a consistant manner
	Input:
		X = x data, can be 1D array or 2D array where each row is a seperate curve
		Y = y data, an array and must have the same dimentions as X
		label = Label for the curves as a string, must have same number of elements as the number of rows of X, (optional)
		xlabel = Label for x-axis as a string
		ylabel = Label for y-axis as a string
		show = boolean, by default False, if set to True, plt.show() will be executed, (optional)
		figure = figure handel to plot to, if None is provieded a new figure will be made, (optional)
		xlim = x-axis limits, [x0, x1], (optional)
		ylim = y-axis limits, [y0, y1], (optional)
		figsize = figure size as a tuple, (optional)
	"""
	if figure == None:
		fig = plt.figure(figsize = figsize);
	else:
		fig = figure;		plt.gcf;

	dim = shape(X);
	if len(dim) == 1:
		plt.plot(X, Y, "x-", label = label);
		plt.grid(True);
		plt.xlim(xlim);
		plt.ylim(ylim);
		plt.xlabel(xlabel);
		plt.ylabel(ylabel);
		plt.legend(loc = "upper left");
	
	elif len(dim) > 1:
		for i in range(dim[0]):
			plt.plot(X[i, :], Y[i, :], "x-", label = label[i]);
		plt.grid(True);
		plt.xlim(xlim);
		plt.ylim(ylim);
		plt.xlabel(xlabel);
		plt.ylabel(ylabel);
		plt.legend(loc = "upper left");

	if show == False:
		return 0;
	elif show == True:
		plt.show();
	return 0;

def Plot_2D(X, Y, Z, xlabel = "You forgot to label your axis", ylabel = "You forgot to label your axis", zlabel = "You forgot to label your axis", figure = None, show = False, xlim = None, ylim = None, zlim = None, figsize = None):
	""" This function will plot Z as a function of X, Y
	Input:
		X = N x N matrix of x nodes (look up documentation of numpy.meshgrid())
		Y = N x N matrix of y nodes (look up documentation of numpy.meshgrid())
		Z = N x N matrix of z values at the x-y nodes
		xlabel = Label for x-axis as a string
		ylabel = Label for y-axis as a string
		zlabel = Label for z-axis as a string
		show = boolean, by default False, if set to True, plt.show() will be executed, (optional)
		figure = figure handel to plot to, if None is provieded a new figure will be made, (optional)
		xlim = x-axis limits, [x0, x1], (optional)
		ylim = y-axis limits, [y0, y1], (optional)
		zlim = z-axis limits, [z0, z1], (optional)
		figsize = figure size as a tuple, (optional)
	"""
	if figure == None:
		fig = plt.figure(figsize = figsize);
	else:
		fig = figure;		plt.gcf;
	ax = plt.axes(projection = "3d");
	ax.plot_surface(X, Y, Z, rstride = 1, cstride = 1, cmap = "jet", edgecolor = "none");
	ax.set_xlabel(xlabel);
	ax.set_ylabel(ylabel);
	ax.set_zlabel(zlabel);

	if show == False:
		return 0;
	elif show == True:
		plt.show();
	return 0;


##%% Example 1 - Simple 1D plot
#N = 50;
#x = linspace(0, 10, N);
#y = random.rand(N);
#label = "Random noise";
#Plot_1D(x, y, show = True, label = label, xlabel = "x [m]");

##%% Example 2 - Multiple curves in one plot
#N = 50;
#x = linspace(0, 10, N);
#M = 3;					# Number of curves
#X = zeros((M, N));
#X[::] = x;
#Y = random.rand(M, N);
#label = array(["Random noise1", "Random noise2", "Random noise3"]);
#Plot_1D(X, Y, label = label, show = True, xlabel = "x [m]", ylabel = "y [m]");
## Change figure size
#size = (10, 5);
#Plot_1D(X, Y, label = label, show = True, xlabel = "x [m]", ylabel = "y [m]", figsize = size);
## Change axis limits
#xlimit = [1, 5];
#Plot_1D(X, Y, label = label, show = True, xlabel = "x [m]", ylabel = "y [m]", figsize = size, xlim = xlimit);

##%% Example 3 - 2D plot
#N = 50;
#x = linspace(-1, 1, N);
#y = linspace(-1, 1, N);
#X, Y = meshgrid(x, y);
#f = lambda x, y: (4*x**3 - 3*x)*(2*y**2 - 1);
#Z = f(X, Y);
#Plot_2D(X, Y, Z, xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]", show = True);



import numpy as np
import matplotlib.pyplot as plt
from Input_parm import *

#Import Airfoils

naca = np.genfromtxt("NASA SC(2)-0606.dat", dtype=str, skip_header=102, skip_footer=0, delimiter="")
naca = naca.astype(np.float)

nacafunc = np.polyfit(naca[:,0],naca[:,1],12)
x = np.linspace(0,1,1001)

def reconstruct_f(coeff, plot_nodes):
	F = 0;
	coeff = coeff[::-1];
	for i in range(len(coeff)):
		F += coeff[i]*plot_nodes**i;
	return F;

y = reconstruct_f(nacafunc,x)

#plt.plot(naca[:,0],naca[:,1])
#plt.plot(x,y)
#plt.show()
naca = np.ones((1001,2))
naca[:,0] = x
naca[:,1] = y

naca2 = np.genfromtxt("NASA SC(2)-0606.dat", dtype=str, skip_header=0, skip_footer=102, delimiter="")
naca2 = naca2.astype(np.float)

nacafunc2 = np.polyfit(naca2[:,0],naca2[:,1],12)
x2 = np.linspace(1,0,1001)

def reconstruct_f(coeff, plot_nodes):
	F = 0;
	coeff = coeff[::-1];
	for i in range(len(coeff)):
		F += coeff[i]*plot_nodes**i;
	return F;

y2 = reconstruct_f(nacafunc2,x2)

plt.plot(naca2[:,0],naca2[:,1])
plt.plot(x2,y2)
plt.show()
naca2 = np.ones((1001,2))
naca2[:,0] = x2
naca2[:,1] = y2

nacatot = np.concatenate((naca2,naca))
np.savetxt('nacatest.dat',nacatot)

print(nacatot)
import numpy as np
import matplotlib.pyplot as plt
from Input_parm import *

#Import Airfoils

naca = np.genfromtxt("NACA 63(2)-615.dat", dtype=str, skip_header=25, skip_footer=0, delimiter="")
naca = naca.astype(np.float)

nasasc = np.genfromtxt("NASA SC(2)-0010.dat", dtype=str, skip_header=102, skip_footer=0, delimiter="")
nasasc = nasasc.astype(np.float)

#Determine chord lengths

mach_area = 0.11965 #from aero textbook, M=0.59
pcw = 0.7 #percentage of span where the strut attaches
y = pcw*span

def Chord_at_y_inb(y,taper_inb,b_inb,S_inb):
	chord_at_y = ((2*S_inb)/((1+taper_inb)*b_inb))*(1-((1-taper_inb)/b_inb)*(y))
	return chord_at_y

if pcw*span <= 36:
  cw = Chord_at_y_inb(y,taper_ratio_inb,span_inboard,S_inb)
  
else:
  print("Fix the chord calculation")

#Verification
y = 36
cwvtip = Chord_at_y_inb(y,taper_ratio_inb,span_inboard,S_inb)
y = 0
cwvroot = Chord_at_y_inb(y,taper_ratio_inb,span_inboard,S_inb)

if (cwvtip - Ct_inb) > 0.01*cwvtip:
	print("Formula Error: tip")

if (cwvroot - Cr_inb) > 0.01*cwvroot:
	print("Formula Error: chord")

else:
	print("good")

cs = 0.45*cw

#translate and scale airoilfs

nasasc = nasasc*cs
naca = naca*cw
nasasc[:,0] = nasasc[:,0]+0.15*cw

#fit curves and make same size lists

nacafunc = np.polyfit(naca[:,0],naca[:,1],12)
x = np.linspace(0,cw,1001)

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
naca [:,1] = -y

nasafunc = np.polyfit(nasasc[:,0],nasasc[:,1],12)
xs = np.linspace(nasasc[0,0],nasasc[-1,0],451)

def reconstruct_f(coeff, plot_nodes):
	F = 0;
	coeff = coeff[::-1];
	for i in range(len(coeff)):
		F += coeff[i]*plot_nodes**i;
	return F;

ys = reconstruct_f(nasafunc,xs)


#plt.plot(naca[:,0],naca[:,1])
#plt.plot(x,y)
#plt.show()
nasasc = np.ones((451,2))
nasasc[:,0] = xs
nasasc[:,1] = -ys


plt.plot(nasasc[:,0],nasasc[:,1])
plt.plot(naca[:,0],naca[:,1])
plt.show()

#Calculate Area between the airfoils

darea = np.zeros(451)
add = np.where(naca[:,0] == cw*0.15)[0][0]
print(nasasc[0])

for i in range(451):
	darea[i] = (naca[(i+add),1] - naca[(add),1])  + nasasc[i,1]

#Calculate combined maximum delta thickness and thus minimum area

biig = darea.max() #maximum delta Area
print("Max delta A=",biig)
indloc = np.where(darea[:] == biig)[0][0]

loc = naca[(add+indloc)]
print("Throat location x/y of Wing (NACA)",loc)
print("Throat location x/y of Strut (NASA SC)",nasasc[indloc])

# Calculate throat where M = 1

throat = biig/mach_area #throat where M=1
print("Throat =", throat)
print(naca[add])

#verification
deltaA = loc[1] - naca[add,1]
VMax_A = deltaA + nasasc[indloc,1]
if (VMax_A - biig) > 0:
	print("Error in area calculations")
else:
	print("good")
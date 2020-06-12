# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 11:52:21 2020

@author: Alberto

This tool combines all calculation tools and plots the von mises stress for
each spanwise station
"""
import numpy as np
from matplotlib import pyplot as plt
from Wingbox_MOI import calc_centroid, calc_MOI
from InterpolationandIntegration import *

# ------------------------  DESIGN INPUTS   ------------------------

#aerodynamics department
wingspan = 50.62/2
c_root = 4.1355
c_tip = 1.81
t_inboard = 0.15 #t/c
t_outboard = 0.18 #t/c
airfoilchange = 18 #m

#material properties
E = 72.4*10**9
G = 28*10**9
v = 0.33 #poisson ratio


#wingbox dimensions
wingbox_start = 0.15
wingbox_end = 0.60
wingbox_height = 0.685 #h/t_airfoil
n_stations = 10

#design parameters
#/!\ make sure each array contains n_station values
# t_top = np.array([0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.002])
# t_bot = np.array([0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.002])
# t_spar = np.array([0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.02,0.002])
# n_stif_top = np.array([5,8,12,15,12,8,5,7,7,2])
# n_stif_bot = np.array([10,16,20,24,15,8,5,5,5,2])
# A_stif = 0.00045
# A_spar_cap = 0.004

t_top = np.array([0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.002])
t_bot = np.array([0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.002,0.002])
t_spar = np.array([0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.02,0.002])
n_stif_top = np.array([4,4,4,4,4,4,4,4,4,4])
n_stif_bot = np.array([4,4,4,4,4,4,4,4,4,4])
A_stif = 0.00045
A_spar_cap = 0.004

#load case
#Mx = np.array([500000,200000,10000])
#Vy = np.array([200000,50000,20000])


# ------------------------  FUNCTIONS   ------------------------

def calc_geometry(wingspan,c_root,c_tip,t_inboard,t_outboard,n_stations):
	"""
	Outputs an array of airfoil chord lengths for n spanwise stations
	Outputs the spanwise y coordinate of n spanwise stations
	"""
	c = np.zeros(n_stations)
	b = np.zeros(n_stations)
	h = np.zeros(n_stations)
	for i in range(n_stations):
		b[i] = i*(wingspan/n_stations)
		c[i] = c_root - (c_root-c_tip)/wingspan * b[i]
		if b[i] < airfoilchange:
			t_c = t_inboard
		else:
			t_c = t_outboard
		h[i] = t_c * c[i]
	return c, b, h

def calc_sigma(Mx,Ixx,z):
	"""
	Outputs an array of normal stresses due to bending
	"""
	sigma = Mx*z/Ixx
	return sigma

def calc_buckling(n_stations, c, n_stif_top, n_stif_bot, t_top, t_bot, sigma_top, sigma_bot, E, v):
	"""
	Outputs the maximum stiffener spacing
	"""
	C = 4 # plates are pinned on 4 sides
	print()
	print("checking for skin buckling...")

	#top skin
	for i in range(n_stations):
		if sigma_top[i] < 0:
			s_top = c[i]/(n_stif_top[i]+1)
			b_max_top = sqrt(C*pi*pi*E*t_top[i]**2/(12*(1-v*v)*abs(sigma_top[i])))
			if s_top > b_max_top:
				print("The top skin buckles at station", i+1,". Increase stiffeners from",n_stif_top[i],"to",np.ceil(c[i]/b_max_top-1))


	#bottom skin
	for j in range(n_stations):
		if sigma_bot[j] < 0:
			s_bot = c[j]/(n_stif_bot[j]+1)
			b_max_bot = sqrt(C*pi*pi*E*t_bot[j]**2/(12*(1-v*v)*abs(sigma_bot[j])))
			if s_bot > b_max_bot:
				print("The bottom skin buckles at station", j+1,". Increase stiffeners from",n_stif_bot[j],"to",np.ceil(c[j]/b_max_bot-1))

	return print("Buckling checked")

def vonmises(sigma,tau):
	"""
	Returns the combined loading stress according to von Misses yeild criterion.
	"""
	A = np.power(sigma,2)
	B = np.power(tau,2)
    #sigm zz and sigma_yy are zero
	return np.sqrt(A/2 + 3*B)

# ------------------------  MAIN   ------------------------

#Compute wingbox geometry
c_airfoil, b, h_airfoil = calc_geometry(wingspan,c_root,c_tip,t_inboard,t_outboard,n_stations)
c = (wingbox_end-wingbox_start) * c_airfoil
h = wingbox_height * h_airfoil

# bending stresses
centroid = calc_centroid(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap)
Ixx, Izz = calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid)


Mx, coeff = RBF_1DInterpol(y_half, M, b)
#print(Mx)


sigma_top = -calc_sigma(Mx,Ixx,h/2-centroid)
sigma_bot = calc_sigma(Mx,Ixx,h/2+centroid)

buckling = calc_buckling(n_stations, c, n_stif_top, n_stif_bot, t_top, t_bot, sigma_top, sigma_bot, E, v)

# verification
print()
print("wingbox widths =",np.round(c,2))
print("spanwise station locations =",b)
print("centroid =", centroid)
print("Ixx, Izz =", Ixx, Izz)
print("Internal Moment =", Mx)

# ------------------------  PLOTS   ------------------------
plt.plot(b, sigma_top, "r--", label = "Normal Stress Top Skin")
plt.plot(b, sigma_bot, "g", label = "Normal Stress Bottom Skin")
plt.ylabel("Stress [N/m^2]")
plt.xlabel("Span [m]")
plt.grid()
plt.yticks(np.arange(-5e8, 5e8, 0.5e8))
plt.legend()
#plt.savefig("Figures/vonmises.png")
plt.show()
#plt.clf()
#plt.close()


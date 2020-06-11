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

# ------------------------  DESIGN INPUTS   ------------------------

#aerodynamics department
wingspan = 50.62/2
c_root = 4.1355
c_tip = 1.81
t_inboard = 0.15 #t/c
t_outboard = 0.18 #t/c
airfoilchange = 18 #m

#wingbox dimensions
wingbox_start = 0.15
wingbox_end = 0.60
wingbox_height = 0.685 #h/t_airfoil
n_stations = 3

#design parameters
# /!\ make sure each array contains n_station values
t_top = np.array([0.004,0.002,0.002])
t_bot = np.array([0.004,0.002,0.002])
t_spar = np.array([0.04,0.02,0.002])
n_stif_top = np.array([20,5,5])
n_stif_bot = np.array([10,5,5])
A_stif = 0.00045
A_spar_cap = 0.004

#load case
Mx = np.array([500000,200000,10000])
Vy = np.array([200000,50000,20000])


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

def calc_buckling(n_stations, c, n_stif_top, n_stif_bot, t_top, t_bot, sigma_top, sigma_bot):
	"""
	Outputs the maximum stiffener spacing
	"""
	C = 4 # plates are pinned on 4 sides
	for i in range(n_stations):
		s_top = c[i]/(n_stif_top[i]+1)
		b_max = sqrt(C*pi*pi*E*t*t/(12*(1-v*v)*sigma))

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

sigma_top = -calc_sigma(Mx,Ixx,h/2-centroid)
sigma_bot = calc_sigma(Mx,Ixx,h/2+centroid)


# verification
print("wingbox widths =",np.round(c,2))
print("spanwise station locations =",b)
print("centroid =", centroid)
print("Ixx, Izz =", Ixx, Izz)


# ------------------------  PLOTS   ------------------------
plt.plot(b, sigma_top, "r--", label = "Compressive stress")
plt.plot(b, sigma_bot, "g", label = "Tensile stress")
plt.ylabel("Stress [N/m^2]")
plt.xlabel("Span [m]")
plt.grid()
plt.yticks(np.arange(-5e8, 5e8, 0.5e8))
plt.legend()
#plt.savefig("Figures/vonmises.png")
plt.show()
#plt.clf()
#plt.close()


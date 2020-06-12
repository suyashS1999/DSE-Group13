# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 11:52:21 2020

@author: Alberto

This tool combines all calculation tools and plots the von mises stress for
each spanwise station
This stress should remain below the material's yield stress
"""
import numpy as np
from matplotlib import pyplot as plt
from Wingbox_MOI import calc_centroid, calc_MOI
from Calculate_Shear_Stresses import Compute_sect_maxShear, V_ , internal_torque
from InterpolationandIntegration import *

# ------------------------  DESIGN INPUTS   ------------------------

#aerodynamics department
wingspan = 50.64/2
c_root = 4.13551818218954
c_tip = 1.81962800016339
t_inboard = 0.14 #t/c
t_outboard = 0.14 #t/c
airfoilchange = 18 #m

#material properties
E = 72.4*10**9
G = 28*10**9
rho_material = 2800 #kg/m^3
v = 0.33 #poisson ratio
tau_allow = 200*10**6


#wingbox dimensions
wingbox_start = 0.15
wingbox_end = 0.60
wingbox_height = 0.685 #h/t_airfoil
n_stations = 10

#design parameters
#/!\ make sure each array contains n_station values

t_top = np.array([0.0005,0.0005,0.001,0.0015,0.002,0.0015,0.0005,0.0025,0.0015,0.0010])
t_bot = np.array([0.001,0.0025,0.0025,0.0025,0.0025,0.0025,0.0015,0.0015,0.0005,0.0005])
t_spar = np.array([0.002,0.003,0.004,0.005,0.005,0.004,0.003,0.002,0.001,0.001])

n_stif_top = np.array([4,4,4,4,4,4,4,21,17,11])
n_stif_bot = np.array([4,16,19,20,19,17,20,10,4,2])

print()
print("Design Inputs are:")
print("stiffener configuration top",n_stif_top)
print("stiffener configuration bot",n_stif_bot)

A_stif = 0.002*0.05*0.04
A_spar_cap = 0.002*0.07*0.05

#load case
#Mx = np.array([500000,200000,10000])
#Vy = np.array([200000,50000,20000])


# ------------------------  FUNCTIONS   ------------------------

def calc_geometry(wingspan,c_root,c_tip,t_inboard,t_outboard,n_stations):
	"""
	Outputs an array of airfoil chord lengths for n spanwise stations
	Outputs the spanwise y coordinate of n spanwise stations
	Outputs the wingbox height
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

def calc_vonmises(sigma,tau_max):
	"""
	Returns the combined loading stress according to von Misses yeild criterion.
	"""
	A = np.power(sigma,2)
	B = np.power(tau_max,2)
    #sigm zz and sigma_yy are zero
	return np.sqrt(A/2 + 3*B)

def calc_mass(c, h, wingspan, n_stations, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap, rho_material):
	Volume = np.zeros(n_stations)
	for i in range(n_stations):
		Area_skin = (t_top[i]+t_bot[i])*c[i]
		Area_spar = t_spar[i]*h[i]
		Area_stif = (n_stif_top[i]+n_stif_bot[i])*A_stif
		Area_spar_caps = 4*A_spar_cap
		Area = Area_skin + Area_spar + Area_stif + Area_spar_caps
		Volume[i] = Area*wingspan/n_stations
	return rho_material*np.sum(Volume)

def calc_fuelvolume(c,h,wingspan,airfoilchange):
	Fuelvol = np.zeros(len(c))
	for i in range(len(c)):
		if i * wingspan/len(c) < airfoilchange:
			Fuelvol[i] = c[i]*h[i]*wingspan/len(c)
	return np.sum(Fuelvol)

# ------------------------  MAIN   ------------------------

#Compute wingbox geometry
c_airfoil, b, h_airfoil = calc_geometry(wingspan,c_root,c_tip,t_inboard,t_outboard,n_stations)
c = (wingbox_end-wingbox_start) * c_airfoil
h = wingbox_height * h_airfoil

# bending stresses
centroid = calc_centroid(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap)
Ixx, Izz = calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid)
Mx, coeff = RBF_1DInterpol(y_half, M, b)
sigma_top = -calc_sigma(Mx,Ixx,h/2-centroid)
sigma_bot = calc_sigma(Mx,Ixx,h/2+centroid)

# shear stress
t_arr = zeros((3,n_stations))
t_arr[0,:] = t_top
t_arr[1,:] = t_bot
t_arr[2,:] = t_spar
tau_max = Compute_sect_maxShear(V, internal_torque, tau_allow, t_inboard, c, t_arr, n_stif_top, n_stif_bot, A_stif, A_spar_cap, iter=False)
print(tau_max)

#vonmises
#vonmises_top = calc_vonmises(sigma_top,tau_max)
#vonmises_bot = calc_vonmises(sigma_bot,tau_max)

#buckling,mass and volume
buckling = calc_buckling(n_stations, c, n_stif_top, n_stif_bot, t_top, t_bot, sigma_top, sigma_bot, E, v)
Mass_wingbox = calc_mass(c, h, wingspan, n_stations, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap, rho_material)
Fuelvolume = calc_fuelvolume(c,h,wingspan,airfoilchange)

#print outputs
print()
#print("wingbox widths =",np.round(c,2))
#print("spanwise station locations =",b)
#print("centroid =", centroid)
#print("Ixx, Izz =", Ixx, Izz)
#print("Internal Moment =", Mx)
print("Normal Stress Top Skin [MPa]",np.round(sigma_top/10**6))
print("Normal Stress Bottom Skin [MPa]",np.round(sigma_bot/10**6))
print("Wingbox mass =",round(Mass_wingbox),"kg (halfwing)")
print("Fuel Volume =",round(Fuelvolume,1),"m^3 (halfwing)")

# ------------------------  PLOTS   ------------------------


ax1 = plt.subplot(221)
#ax1.set_title("Thickness Distribution")
ax1.set_ylabel("Thickness [mm]")
ax1.set_xlabel("Span [m]")
ax1.step(b, t_top*1000, "r", label = "Top skin")
ax1.step(b, t_bot*1000, "b", label = "Bottom skin")
ax1.step(b, t_spar*1000, "g", label = "Spar")
ax1.grid()
ax1.set_yticks(arange(0, 6, 1))
ax1.tick_params(axis='y')
ax1.legend(loc='upper left')

ax2 = plt.subplot(222)
#ax2.set_title("Stiffener Distribution")
ax2.set_ylabel("Number of stiffeners")
ax2.set_xlabel("Span [m]")
ax2.step(b, n_stif_top, "r", label = "Top skin")
ax2.step(b, n_stif_bot, "b", label = "Bottom skin")
ax2.grid()
ax2.set_yticks(arange(0, 30, 5))
ax2.tick_params(axis='y')
ax2.legend(loc='upper left')

ax3 = plt.subplot(212)
#ax3.set_title("Stress Distribution")
ax3.set_ylabel("Stress [MPa]")
ax3.set_xlabel("Span [m]")
ax3.plot(b, sigma_top/10**6, "r", label = "Top Skin")
ax3.plot(b, sigma_bot/10**6, "b", label = "Bottom Skin")
ax3.grid()
ax3.set_yticks(arange(-500, 500, 100))
ax3.tick_params(axis='y')
ax3.legend(loc='upper left')


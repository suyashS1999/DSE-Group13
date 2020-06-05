# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:13:41 2020

@author: Alberto

This script computes the Moment of Inertia along n stations of a simple
rectangular wingbox using two spars, top and bottom skin, four spar caps
and n stiffeners
"""

import numpy as np
from math import *


def calc_centroid(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap):
	"""Return the z coordinate of the centroid. (z points upward)"""

	sumAreaZ = n_stif_top*A_stif*h/2 + n_stif_bot*A_stif*(-h/2)
	sumArea = n_stif_top*A_stif + n_stif_bot*A_stif+ 2*t_spar*h + t_top*c + t_bot*c
	z_centroid = sumAreaZ/sumArea

	return z_centroid


def calc_stif_locations(c, h, n_stif_top, n_stif_bot):
	""" Function to compute coordinates of all the top and bottom stiffeners
		Input:
			c: width
			h: heighth
		Output:
			stif_coordinates: 2D array with first row having the top stiffener coordinates with respect
			to wingbox centre and second row having the bottom stiffener coordinates
			[-c,-c+0.005,...]
	"""
	stif_coordinates_top = np.ones((len(n_stif_top),len(c)))
	stif_coordinates_bot = np.ones((len(n_stif_bot),len(c)))

	#s_top = stiffener spacing on top skin
	s_top = np.ones(len(c))
	s_bot = np.ones(len(c))
	for i in range(len(c)):
		s_top[i] = c[i]/(n_stif_top[i]+1) #remember there are 2 spar caps in the corners, so
		s_bot[i] = c[i]/(n_stif_bot[i]+1) #if you have 2 stiffeners there are 3 spacings

	for i in range(len(n_stif_top)):
		stif_coordinates_top[:,i] = -c/2 + s_top[i]*(i+1) #index plus one because start with spar cap

	for j in range(len(n_stif_bot)):
		stif_coordinates_bot[:,j] = -c/2 + s_bot[j]*(j+1) #index plus one because start with spar cap

	return stif_coordinates_top, stif_coordinates_bot

def calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid,stif_coordinates_top,stif_coordinates_bot):
	""" Function to compute MOI of rectangular wingbox
		Input:
			c: width
			h: heighth
		Output:
			MOI: array of MOI for n stations
	"""
	MOI = np.ones(len(t_top))

	# ------------------------  I_xx  ------------------------
	# I_xx consists of 4 parts: spars (1), skin plates (2), stiffeners (3), spar caps (4)

	I_xx1 = 2 * 1/12*t_spar*h**3
	I_xx2 = 1/12*c*t_top**3 + c*t_top*(h/2-centroid)**2 + 1/12*c*t_bot**3 + c*t_bot*(-h/2-centroid)**2
	I_xx3 = n_stif_top * A_stif*(h/2-centroid)**2 + n_stif_bot * A_stif*(-h/2-centroid)**2
	I_xx4 = 2 * A_spar_cap*(h/2-centroid)**2 + 2 * A_spar_cap*(-h/2-centroid)**2

	Ixx = I_xx1 + I_xx2 + I_xx3 + I_xx4

	# ------------------------  I_zz  -------------------------
	# I_yy consists of 4 parts: spars (1), skin plates (2), stiffeners (3), spar caps (4)

	I_zz1 = 2 * ( 1/12*h*t_spar**3 + t_spar*(c/2)**2 )
	I_zz2 = 1/12*t_top*c**3 + 1/12*t_bot*c**3
	I_zz3 = 0
	I_zz4 = 2 * A_spar_cap*(c/2)**2 + 2 * A_spar_cap*(-c/2)**2

	Izz = I_zz1 + I_zz2 + I_zz3 + I_zz4

	return Ixx, Izz



# ------------------------  TEST CASE   ------------------------

c = np.array([4,2,2])
h = np.array([4,2,2])
t_top = np.array([0.004,0.002,0.002])
t_bot = np.array([0.004,0.002,0.002])
t_spar = np.array([0.04,0.02,0.002])
n_stif_top = np.array([20,5,5])
n_stif_bot = np.array([10,5,5])
A_stif = 0.00045
A_spar_cap = 0.004


centroid = calc_centroid(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap)
stif_coordinates_top, stif_coordinates_bot = calc_stif_locations(c, h, n_stif_top, n_stif_bot)
Ixx, Izz = calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid,stif_coordinates_top,stif_coordinates_bot)

print("centroid equals", centroid)
print("stiffener coordinates on skin",stif_coordinates_top, stif_coordinates_bot)
print("Ixx, Izz", Ixx, Izz)
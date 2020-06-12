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

	sumAreaZ = n_stif_top*A_stif*h/2 + n_stif_bot*A_stif*(-h/2) + c*t_top*h/2 + c*t_bot*(-h/2)
	sumArea = n_stif_top*A_stif + n_stif_bot*A_stif+ 2*t_spar*h + t_top*c + t_bot*c
	z_centroid = sumAreaZ/sumArea

	return z_centroid


def calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid):
	""" Function to compute MOI of rectangular wingbox
		Input:
			c: width
			h: heighth
		Output:
			MOI: array of MOI for n stations
	"""

	# ------------------------  I_xx  ------------------------
	# I_xx consists of 4 parts: spars (1), skin plates (2), stiffeners (3), spar caps (4)

	I_xx1 = 2 * 1/12*t_spar*h**3
	#print (I_xx1)
	I_xx2 = 1/12*c*t_top**3 + c*t_top*(h/2-centroid)**2 + 1/12*c*t_bot**3 + c*t_bot*(-h/2-centroid)**2
	#print (I_xx2)
	I_xx3 = n_stif_top * A_stif*(h/2-centroid)**2 + n_stif_bot * A_stif*(-h/2-centroid)**2
	#print (I_xx3)
	I_xx4 = 2 * A_spar_cap*(h/2-centroid)**2 + 2 * A_spar_cap*(-h/2-centroid)**2
	#print (I_xx4)

	Ixx = I_xx1 + I_xx2 + I_xx3 + I_xx4

	# ------------------------  I_zz  -------------------------
	# I_yy consists of 4 parts: spars (1), skin plates (2), stiffeners top (3), stiffeners bottom (4), spar caps (5)

	I_zz1 = 2 * ( 1/12*h*t_spar**3 + t_spar*(c/2)**2 )
	I_zz2 = 1/12*t_top*c**3 + 1/12*t_bot*c**3

	# to find I_zz3, the stiffener locations must be known
	MOI = np.zeros(len(n_stif_top))
	for i in range(len(n_stif_top)):
		A_stif_top = A_stif * np.ones(n_stif_top[i])
		s_top = c[i]/(n_stif_top[i]+1) #remember there are 2 spar caps in the corners, so
		stif_coordinates_top = np.zeros(n_stif_top[i])
		#print(stif_coordinates_top)
		for k in range(n_stif_top[i]):
			stif_coordinates_top[k] = -c[i]/2 + s_top*(k+1) #index plus one because start with spar cap
		MOI[i] = A_stif_top.dot(stif_coordinates_top**2)
	I_zz3 = MOI
	#print(I_zz3)

	# I_zz4 (stiffners bottom)
	MOI = np.zeros(len(n_stif_bot))
	for i in range(len(n_stif_bot)):
		A_stif_bot = A_stif * np.ones(n_stif_bot[i])
		s_bot = c[i]/(n_stif_bot[i]+1) #remember there are 2 spar caps in the corners, so
		stif_coordinates_bot = np.zeros(n_stif_bot[i])
		for k in range(n_stif_bot[i]):
			stif_coordinates_bot[k] = -c[i]/2 + s_bot*(k+1) #index plus one because start with spar cap
		MOI[i] = A_stif_bot.dot(stif_coordinates_bot**2)
	I_zz4 = MOI
	#print(I_zz4)

	I_zz5 = 2 * A_spar_cap*(c/2)**2 + 2 * A_spar_cap*(-c/2)**2

	Izz = I_zz1 + I_zz2 + I_zz3 + I_zz4 + I_zz5

	return Ixx, Izz



# ------------------------  TEST CASE   ------------------------

# c = np.array([4,2,2])
# h = np.array([4,2,2])
# t_top = np.array([0.004,0.002,0.002])
# t_bot = np.array([0.004,0.002,0.002])
# t_spar = np.array([0.04,0.02,0.002])
# n_stif_top = np.array([20,5,5])
# n_stif_bot = np.array([10,5,5])
# A_stif = 0.00045
# A_spar_cap = 0.004


# centroid = calc_centroid(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap)
# Ixx, Izz = calc_MOI(c, h, t_top, t_bot, t_spar, n_stif_top, n_stif_bot, A_stif, A_spar_cap,centroid)

# print("centroid equals", centroid)
# print("Ixx, Izz", Ixx, Izz)

# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 10:22:17 2020

@author: Alberto

This script computes the loading on the hinge of a folding wingtip.
The loading is found for three cases:
1) ground gusts when folded
2) maximum loadfactor at manoevre speed
3) Maybe a third case, a birdstrike of a 4lb bird (if I get to that).

This script enables finding the motor torque as well as enables sizing the lugs and locking pins made out of steel.
"""

import numpy as np
from math import *


# ----------------------- AERODYNAMIC SPANWISE LIFT -------------------------
# L_prime is the Cl multiplied with c/cref already
L_prime = np.array([0.45401197, 0.59314027, 0.67983843, 0.73606313, 0.77251623, 0.80735009,
 0.8353416,  0.86410226, 0.8862658,  0.91311351, 0.93534546, 0.96534544,
 1.01009516, 0.93205674, 0.99882457, 1.05348972, 1.10355137, 1.14644941,
 1.18648774, 1.22565622, 1.24058785, 1.26768158, 1.27972376, 1.28495111,
 1.27576291, 1.22487137, 1.10590927, 1.10590927, 1.22487137, 1.27576291,
 1.28495111, 1.27972376, 1.26768158, 1.24058785, 1.22565622, 1.18648774,
 1.14644941, 1.10355137, 1.05348972, 0.99882457, 0.93205674, 1.01009516,
 0.96534544, 0.93534546, 0.91311351, 0.8862658 , 0.86410226, 0.8353416,
 0.80735009, 0.77251623, 0.73606313, 0.67983843, 0.59314027, 0.45401197])

spanlocation = np.array([-24.35585, -23.84732, -23.33879, -22.83026, -22.32173, -21.8132, -21.30466,
 -20.79613, -20.28759, -19.77905, -19.27052, -18.76198, -18.25344, -17.35214,
 -16.06664, -14.78112, -13.49558, -12.21003, -10.92447,  -9.6389 ,  -8.35332,
  -7.06773,  -5.78213,  -4.49653,  -3.21092,  -1.9253 ,  -0.63968,   0.63968,
   1.9253,    3.21092,   4.49653,   5.78213,   7.06773,   8.35332,   9.6389,
  10.92447,  12.21003,  13.49558,  14.78112,  16.06664,  17.35214,  18.25344,
  18.76198,  19.27052,  19.77905,  20.28759,  20.79613,  21.30466,  21.8132,
  22.32173,  22.83026,  23.33879,  23.84732,  24.35585])

# ----------------------------- FUNCTIONS -----------------------------------
def ISA_trop(h):
	""" This function computes the atmospheric properties
		within the troposphere
	Input:
		h = altitude [m]
	Output:
		T = Temperature [K]
		p = Pressure [Pa]
		rho = Densituy [kg/m^3]
		a = spped of sound [m/s]
	"""
	T = 288.15 - 0.0065*h;
	p = 101325*(T/288.15)**(-g/(-0.0065*287));
	rho = 1.225*(T/288.15)**(-g/(-0.0065*287) - 1);
	a = sqrt(1.4*287*T);
	return T, p, rho, a;

def stationspan(spanlocation,b):
	# Verified by Alberto with sample 2,4,4,2 chord distribution
	""" This function computes surface area per spanwise station,
	where Cl is applied in the middle of that section
	Input:
		spanlocation: spanwise location of each station [m]
		b: wingspan [m]
	Output:
		sectionsurface: array of surface areas where each section CL applies to
	"""
	n_stations = len(chordlist)
	b_prime = b/n_stations
	Sectionsurface = chordlist * b_prime
	return Sectionsurface


# ------------------------------ INPUTS -------------------------------------
g = 9.81
h = 11000
v_manoeuvre = 100
S = 142.5
MAC = 1.5
b = 2*24.61
n_ult = 3.75
T, p, rho_cruise, a = ISA_trop(h)
print(rho_cruise)


# Compute Lift from CL
L = np.ones((len(spanlocation)))

b_prime = b/len(spanlocation)
for i in range(len(spanlocation)):
	L[i] = L_prime[i]*0.5*1.225*v_manoeuvre**2*b_prime*MAC

print(sum(L))



# Load case: load factor of 3.75 at manoevre speed
L[i] = 3.75*L[i]

# The wingtip runs from 36.0m to 24.61m
# Only the aerodynamic loads there are needed, so others are filtered out




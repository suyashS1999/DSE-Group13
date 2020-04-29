# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from numpy import*

def RangeJet(V, cj, L_D, Weight_frac):
	""" This function computes the range that can be obtained by a 
		jet aircraft
	Input:
		V = cruise velocity [m/s]
		cj = Thrust specific fuel consumption [kg/Ns]
		L_D = Lift to Drag ratio
		Weight_frac = weight fraction before and after cruise phase [-]
	Output:
		R = Range [m]
	"""
	g = 9.80665;
	R = (V/(cj*g))*L_D*ln(Weight_frac);
	return R;

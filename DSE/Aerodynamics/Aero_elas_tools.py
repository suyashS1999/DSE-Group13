# -*- coding: utf-8 -*-
"""

"""

import numpy as np
from matplotlib import pyplot as plt




def diverg_q(K_th,K_phi,S,e,a,sweep_half,b):
	"""Semi rigid wing model (including sweep)[from Weisshaar pg 110 of book]"""
	
	v1 = (K_th/S*e*a)
	v2 = (np.cos(sweep_half)**2)*(1-((b/e)*(K_th/K_phi)*0.5*np.tan(sweep_half)))
	
	q = v1/v2 
	
	return q


def crit_sweep(e,c,b,K_th,K_phi):
	
	lmb = np.arctan(2*(e/c)*(c/b)*(K_phi/K_th))
	
	return np.rad2deg(lmb)


def aileron_eff_full(rho,V_inf,a,q_d,K_th,cm_d,cl_d,c,S,sweep_half):
	
	"Assuming full span aileron [from weisshaar pg 117 ]"
	
	q = 0.5*rho*(V_inf**2)
	
	v1 = 1+(((q*S*a*c*cm_d)/(K_th*cl_d))*(np.cos(sweep_half)**2))
	v2 = 1-(q/q_d)
	
	return v1/v2

def aileron_eff_act(lamb,half_span,y_a,ec,cl_a,cl_d,c,cm_d):
	
	"""Using semi rigid wing model from book
	lamb - see below fucntion
	half_span
	y_a - postion of aileron (starting) on half span
	ec - see figure
	cl_a - 2D lift coeff
	cl_d - 2D aileron life coeff w.r.t aileron def
	c - mean chord
	cm_d - 2D change of moment wrt delta"""
	
	e = ec/c
	
	v1 = ((np.cos(lamb*y_a)/np.cos(lamb*half_span))-1)*(1/cl_a)*(cl_d)
	
	v2 = ((np.cos(lamb*y_a)/np.cos(lamb*half_span))-1-(((lamb**2)*((half_span**2)-(y_a**2))*0.5)))
	
	v3 = (1/(e*cl_a))*cm_d
	
	v4 = ((np.tan(lamb*half_span)/(lamb*half_span)) - 1)
	
	eff = (v1 + (v2*v3))/v4
	
	return eff
	
def lamb(rho,V,ec,cl_a,G,J):
	"""rho - density
	V - free streamm vel
	ec - distance between aero center and flexural axis
	cl_a - 2D lift coeff
	G - 
	J - """
	
	
	res = np.sqrt((0.5*rho*(V**2)*(ec**2)*cl_a)/(G*J))
	
	return res


	

def diverg_speed_straight(CL_alp,G,J,rho,ec,half_span):
	
	"""Diverg speed for straight wing (no sweep) from book"""
	
	vel = np.sqrt(((np.pi**2)*G*J)/(2*rho*(ec**2)*(half_span**2)*CL_alp))
	
	return vel


			   
   




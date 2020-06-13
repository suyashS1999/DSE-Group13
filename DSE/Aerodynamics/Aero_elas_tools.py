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


def cl_d(a0,E):
	
	val = (a0/np.pi)*((np.arccos(1-(2*E))) + (2*np.sqrt(E*(1-E))))	
	
	return val

def cm_d(a0,E):
	
	val = -(a0/np.pi)*(1-E)*(np.sqrt(E*(1-E)))
	
	return val

def diverg_speed_straight(CL_alp,G,J,rho,ec,half_span):
	
	"""Diverg speed for straight wing (no sweep) from book"""
	
	vel = np.sqrt(((np.pi**2)*G*J)/(2*rho*(ec**2)*(half_span**2)*CL_alp))
	
	return vel


def A_matrix(m,half_span,c,x_f):
	
	A = np.zeros((2,2))
	
	A[0,0] = half_span*c/5
	A[0,1] = (half_span/4)*(((c**2)*0.5)-(c*x_f))
	A[1,0] = (half_span/4)*(((c**2)*0.5)-(c*x_f))
	A[1,1] = (half_span/3)*(((c**3)/3) - ((c**2)*x_f) + (c*(x_f**2)))
	
	return A*m

def B_matrix(c,half_span,aw,e,M_theta_dot):
	
	B = np.zeros((2,2))
	B[0,0] = (c*half_span*aw)/10
	B[1,0] = -((c**2)*half_span*e*aw)/8
	B[1,1] = -((c**3)*half_span*M_theta_dot)/24

	return B

def C_matrix(c,half_span,aw,e):
	C = np.zeros((2,2))	
	
	C[0,1] = (c*half_span*aw)/8
	C[1,1] = -((c**2)*half_span*e*aw)/6
	
	return C

def E_matrix(EI,GJ,half_span):
	E= np.zeros((2,2))
	
	E[0,0] = (4*EI)/(half_span**3)
	E[1,1] = (GJ)/half_span
	
	return E


	
		   
   



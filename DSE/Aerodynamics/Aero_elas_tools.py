# -*- coding: utf-8 -*-
"""

"""

import numpy as np
from matplotlib import pyplot as plt




def diverg_q(K_th,K_phi,half_area,e,a,c,sweep_half,b):
	"""Semi rigid wing model (including sweep)[from Weisshaar pg 110 of book]"""
	
	v1 = (K_th/half_area*e*c*a)
	v2 = (np.cos(sweep_half)**2)*(1-((b/e*c)*(K_th/K_phi)*0.5*np.tan(sweep_half)))
	
	q = v1/v2 
	
	return q


def diverg_q_v2(GJ,EI,l,e,a,c,sweep_half):
 	"""Semi rigid wing model (including sweep)[from New book Into to str dynamic and aeroelas]"""
 	
 	v1 = (GJ*(np.pi**2))
 	v2 = (4*e*c*c*a*(l**2)*(np.cos(sweep_half)**2))*(1-(((3*(np.pi**2)*l)/(76*e*c))*(GJ/EI)*np.tan(sweep_half)))
 	
 	q = v1/v2 
 	
 	return q


def crit_sweep(EI,GJ,e,c,l):
	
	lmb = np.arctan((76*EI*e*c)/(3*(np.pi**2)*GJ*l))
	
	return np.rad2deg(lmb)


def aileron_eff_full(rho,V_inf,a,q_d,K_th,cm_d,cl_d,c,S,sweep_half):
	
	"Assuming full span aileron [from weisshaar pg 117 ]"
	
	q = 0.5*rho*(V_inf**2)
	
	v1 = 1+(((q*S*a*c*cm_d)/(K_th*cl_d))*(np.cos(sweep_half)**2))
	v2 = 1-(q/q_d)
	
	return v1/v2

def aileron_eff_act(lamb,half_span,y_a,e,cl_a,cl_d,c,cm_d):
	
	"""Using semi rigid wing model from book
	lamb - see below fucntion
	half_span
	y_a - postion of aileron (starting) on half span
	e - see figure
	cl_a - 2D lift coeff
	cl_d - 2D aileron life coeff w.r.t aileron def
	c - mean chord
	cm_d - 2D change of moment wrt delta"""
	

	
	v1 = ((np.cos(lamb*y_a)/np.cos(lamb*half_span))-1)*(1/cl_a)*(cl_d)
	
	v2 = ((np.cos(lamb*y_a)/np.cos(lamb*half_span))-1-(((lamb**2)*((half_span**2)-(y_a**2))*0.5)))
	
	v3 = (1/(e*cl_a))*cm_d
	
	v4 = ((np.tan(lamb*half_span)/(lamb*half_span)) - 1)
	
	eff = (v1 + (v2*v3))/v4
	
	return eff
	
def lamb(rho,V,e,c,cl_a,GJ):
	"""rho - density
	V - free streamm vel
	ec - distance between aero center and flexural axis
	cl_a - 2D lift coeff
	G - 
	J - """
	
	
	res = np.sqrt((0.5*rho*(V**2)*(e*c**2)*cl_a)/(GJ))
	
	return res


def cl_d(a0,E):
	
	val = (a0/np.pi)*((np.arccos(1-(2*E))) + (2*np.sqrt(E*(1-E))))	
	
	return val

def cm_d(a0,E):
	
	val = -(a0/np.pi)*(1-E)*(np.sqrt(E*(1-E)))
	
	return val

def diverg_speed_straight(a,GJ,rho,e,c,half_span):
	
	"""Diverg speed for straight wing (no sweep) from book"""
	
	q = ((np.pi**2)*GJ)/(4*(e*c**2)*(half_span**2)*a)
	
	return q


def aileron_eff_new(lamb,l,c,cm_d,cl_d,e,a):
	
	v1 = lamb*l*((c*cm_d*(((lamb*l)**2) - (2*(1/np.cos(lamb*l))) +2)) - (2*e*c*cl_d*((1/np.cos(lamb*l)) -1)))
	v2 = 2*a*e*c*((lamb*l) - np.tan(lamb*l))
	
	return v1/v2
	
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


	
		   
   



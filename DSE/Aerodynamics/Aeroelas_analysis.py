# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import Aero_elas_tools as tool
from matplotlib import pyplot as plt
# import Input_parm as inp
# from Aero_tools import *




c = 2
# c = ave_chord      # Mean chord  
# e = 0.25*ave_chord # distance between ac line and elastic line
# half_span = span/2  # half span
half_span = 7.5
# sweep_half = sweep_x(0.5,sweep_LE,C_r_m,C_t_m,span) # radians
y_a = 10.8
x_f = 0.48*c
e = (x_f/c) - (1/4)
m = 200# mass/area
M_theta_dot = -1.2

K_th = 1# Torsional stiffness
K_phi = 1# Bending stiffness

# CL_alp = CL_alpha_DATCOM(AR,M_cruise,a0_full,sweep_half) #wing life slope
GJ = 2e6
EI = 2e7


E_f = 0.4  #control_chord/total_chord
# rho = rho_cruise
rho = 1.225
# cl_a = a0_full
cl_a = 2*np.pi
# cm_d = tool.cl_d(cl_a,E_f)
# cl_d = tool.cm_d(cl_a,E_f)



# v_div_str = tool.diverg_speed_straight(CL_alp,G,J,rho,e,half_span)

# v_div_swept = tool.diverg_q(K_th,K_phi,S,e,CL_alp,sweep_half,half_span)

# crit_sweep = tool.crit_sweep(e,c,half_span,K_th,K_phi)

V_inf = np.linspace(10,170,1000)

# q_d = 0.5*rho*(v_div_swept**2)
# ail_eff_full_s = tool.aileron_eff_full(rho,V_inf,CL_alp,q_d,K_th,cm_d,cl_d,c,S,sweep_half)

# lamb = tool.lamb(rho,V_inf,e,cl_a,G,J)

# ail_eff_str_wing = tool.aileron_eff_act(lamb,half_span,y_a,e,cl_a,cl_d,c,cm_d)


dummy_1 = np.zeros((2,2))
I_dummy = np.eye(2)

A_mat = tool.A_matrix(m, half_span, c, x_f)
# B_mat = tool.B_matrix(c, half_span, cl_a, e, M_theta_dot)
B_mat = np.zeros((2,2))
C_mat = tool.C_matrix(c, half_span, cl_a,e)
E_mat = tool.E_matrix(EI, GJ, half_span)

freq = np.zeros((len(V_inf),4))
damp = np.zeros((len(V_inf),4))

freq_2 = np.zeros((len(V_inf),2))
damp_2 = np.zeros((len(V_inf),2))


for i in range(len(V_inf)):
	V = V_inf[i]
	A_1 = (rho*(V**2)*C_mat)+E_mat
	A_2 = rho*V*B_mat
	Q = np.block([[dummy_1,I_dummy],[-np.linalg.inv(A_mat).dot(A_1),-np.linalg.inv(A_mat).dot(A_2)]])
	check= np.linalg.eigvals(Q)

	
	for j in range(4):
		a = check[j].real
		b = check[j].imag
		
		
		freq[i,j] = (np.sqrt((a**2)+(b**2)))
		damp[i,j] = -100*a/(freq[i,j])
		freq[i,j] = (np.sqrt((a**2)+(b**2)))/(2*np.pi)
	



print("Results \n \n")

# print("Diverg Speed (Swept wing model)", v_div_swept,"\n")
# print("Diverg Speed (Straight Wing", v_div_str,"\n")
# print("Critical Sweep", crit_sweep,"\n")



# plt.figure()
# plt.plot(V_inf,ail_eff_full_s,label="full span aileron")
# plt.plot(V_inf,ail_eff_str_wing,label="Actual aileron (no sweep and 2D coeff)")
# plt.xlabel("V_inf")
# plt.ylabel("Äileron effectiveness")
# plt.grid()

plt.figure()
plt.subplot(2,1,1)
plt.plot(V_inf,damp,label="Okeii")
plt.xlabel("V_inf")
plt.ylabel("Damping Ratio [%]")
plt.legend()
plt.subplot(2,1,2)
plt.plot(V_inf,freq)
plt.xlabel("V_inf")
plt.ylabel("Frequency")
plt.title("Negative")

# plt.figure()
# plt.subplot(2,1,1)
# plt.plot(V_inf,damp_2[:,0])
# plt.plot(V_inf,damp_2[:,1])
# plt.xlabel("V_inf")
# plt.ylabel("Damping Ratio [%]")

# plt.subplot(2,1,2)
# plt.plot(V_inf,freq_2[:,0])
# plt.plot(V_inf,freq_2[:,1])
# plt.xlabel("V_inf")
# plt.ylabel("Frequency")
# plt.title("Postitive")
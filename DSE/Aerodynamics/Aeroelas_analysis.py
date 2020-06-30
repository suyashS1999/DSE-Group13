# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import Aero_elas_tools as tool
from matplotlib import pyplot as plt
import Input_parm as inp
from Aero_tools import *





c =  inp.MAC[0]    # Mean chord  
e = 0.22 # distance between ac line and elastic line
half_span = inp.span/2  # half span
half_area = inp.S*0.5
sweep_half = sweep_x(0.5,inp.sweep_LE,inp.C_r_m,inp.C_t_m,inp.span) # radians
y_a = 15.5
x_f = 0.5*c
m = 62.1733 # mass/area
M_theta_dot = -1.2



CL_alp = CL_alpha_DATCOM(inp.AR,inp.M_cruise,inp.a0_full,sweep_half) #wing life slope
GJ = 28e9 *0.001366181409646969 #  0.0007449479615522029 #0.001166181409646969 # # 
EI = 72e9 * 0.00029890223261558427 
 

E_f = 0.4  #control_chord/total_chord
rho = inp.rho_cruise
cl_a = inp.a0_full
cm_d = tool.cm_d(cl_a,E_f)
cl_d = tool.cl_d(cl_a,E_f)
# CL_alp = cl_a


v_div_str = (tool.diverg_speed_straight(CL_alp,GJ,rho,e,c,half_span))

q_d = tool.diverg_q_v2(GJ,EI,half_span,e,cl_a,c,sweep_half)
# 
crit_sweep = tool.crit_sweep(EI,GJ,e,c,half_span)

verifi = v_div_str - tool.diverg_q_v2(GJ,EI,half_span,e,CL_alp,c,0)

V_inf = np.linspace(0,300,100)

# ail_eff_full_s = tool.aileron_eff_full(rho,V_inf,CL_alp,q_d,K_th,cm_d,cl_d,c,half_area,sweep_half)

lamb = tool.lamb(rho,V_inf,e,c,cl_a,GJ)

ail_eff_str_wing = tool.aileron_eff_act(lamb,half_span,y_a,e,cl_a,cl_d,c,cm_d)

ail_eff_v2 = tool.aileron_eff_new(lamb,half_span,c,cm_d,cl_d,e,cl_a)

dummy_1 = np.zeros((2,2))
I_dummy = np.eye(2)

A_mat = tool.A_matrix(m, half_span, c, x_f)
B_mat = tool.B_matrix(c, half_span, cl_a, e, M_theta_dot)
# B_mat = np.zeros((2,2))
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

print("Diverg Speed (q_d) (Swept wing model)", q_d,"\n")
print("Diverg Speed (q_d) (Straight Wing)", v_div_str,"\n")
print("Critical Sweep", crit_sweep,"\n")
print("Verification =",verifi)

plt.rcParams.update({'font.size': 30})

plt.figure(num=None, figsize=(15, 12), dpi=80)
# plt.plot(V_inf,ail_eff_v2,label="Actual aileron (no sweep and 2D coeff)")
plt.plot(V_inf,ail_eff_str_wing,label="Actual aileron (v2)",color="red",linewidth=2.0)
plt.xlabel("$V_{\infty}$ [m/s]")
plt.ylabel("Aileron effectiveness [" +"$pb_{1/2}/V_{\infty} \beta $" +"]")
# plt.ylim((0.7,1.1))
# plt.legend()
plt.grid()
plt.show()

plt.figure(num=None, figsize=(12, 11), dpi=80)
# plt.figure()
plt.plot(V_inf,damp[:,1],label="Bending Mode",color="red",linewidth=2.0)
plt.plot(V_inf,damp[:,2],label="Torsion Mode",color="blue",linewidth=2.0)
plt.axvline(290,linestyle="dashed", color='black',label="Critical Flutter Speed",linewidth=2.0)
plt.xlabel("$V_{\infty} [m/s]$")
plt.ylabel("Damping Ratio [%]")
plt.ylim((-5,105))
plt.legend()
plt.grid()

plt.figure(num=None, figsize=(12, 11), dpi=80)
# plt.figure()
plt.plot(V_inf,freq[:,1],label="Bending Mode",color="red",linewidth=2.0)
plt.plot(V_inf,freq[:,2], label="Torsion Mode",color="blue",linewidth=2.0)
plt.axvline(290,linestyle="dashed",color='black', label="Critical Flutter Speed",linewidth=2.0)
plt.xlabel("$V_{\infty} [m/s]$")
plt.ylabel("Frequency [Hz]")
plt.ylim((0,23))
plt.grid()
plt.legend(loc="upper left")
plt.show()

# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import init_tools as tool
from matplotlib import pyplot as plt


S = 100  # Wing surface area [m^2]
e = 0.4  # distance between ac line and elastic line
c = 4.2  # Mean chord  
half_span = 25   # half span
sweep_half = 0.5 # radians
y_a = 10

K_th = 1# Torsional stiffness
K_phi = 1# Bending stiffness

CL_alp = 5 #wing life slope
G = 1
J = 1
rho = 1.225
cm_d = 1
cl_d = 1
cl_a = 1


v_div_str = tool.diverg_speed_straight(CL_alp,G,J,rho,e,half_span)

v_div_swept = tool.diverg_q(K_th,K_phi,S,e,CL_alp,sweep_half,half_span)

crit_sweep = tool.crit_sweep(e,c,half_span,K_th,K_phi)

V_inf = np.linspace(1,250,100)

q_d = 0.5*rho*(v_div_swept**2)
ail_eff_full_s = tool.aileron_eff_full(rho,V_inf,CL_alp,q_d,K_th,cm_d,cl_d,c,S,sweep_half)

lamb = tool.lamb(rho,V_inf,e,cl_a,G,J)

ail_eff_str_wing = tool.aileron_eff_act(lamb,half_span,y_a,e,cl_a,cl_d,c,cm_d)



print("Results \n \n")

print("Diverg Speed (Swept wing model)", v_div_swept,"\n")
print("Diverg Speed (Straight Wing", v_div_str,"\n")
print("Critical Sweep", crit_sweep,"\n")



plt.figure()
plt.plot(V_inf,ail_eff_full_s,label="full span aileron")
plt.plot(V_inf,ail_eff_str_wing,label="Actual aileron (no sweep and 2D coeff)")
plt.xlabel("V_inf")
plt.ylabel("Äileron effectiveness")
plt.grid()
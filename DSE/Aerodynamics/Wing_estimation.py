import numpy as np
import Aero_tools as tl
from Input_parm import *
import time as time

#%% ------------------ Main ----------------------
t_s = time.time();

sweep_half_c = tl.sweep_x(0.5, sweep_LE, C_r_m, C_t_m, span);			# Half chord sweep
sweep_quart = tl.sweep_x(0.25, sweep_LE, C_r_m, C_t_m, span);			# Quarter chord sweep
CLmax_Clmax = 0.73;

Alpha0_inb = -4;					# Zero lift aoa inboard airfoil [deg]
Alpha0_outb = -2;					# Zero lifr aoa outboard airfoil [deg]
d_alpha_clmax = 3.5					# [deg]

# ----- Reading graphs from NACA tests ------

V_M02 = 59.01;							# [m/s]
kin_visc = 1.46e-05						# kinematic viscoty at sea level
MAC = MAC[0];							# Mean aerodynamic chord [m]	
ReM02 = (V_M02*MAC)/kin_visc			# Reynolds number at MAC

if ReM02 >= 8.1e06:
	Clmax_M02_inb = 1.67;
	Clmax_M02_outb = 1.6;

elif ReM02 >= 5.58e06 and ReM02 < 8.1e06:
	Clmax_M02_inb = 1.6;
	Clmax_M02_outb = 1.5;

elif ReM02 < 5.58e06:
	Clmax_M02_inb = 1.47;
	Clmax_M02_outb = 1.38;
else:
	print("Clmax at M = 0.2 is wronk");


# ----- Lift Calculations -------
# INBOARD
CL_alp_inb = tl.CL_alpha_DATCOM(AR_inb, M_cruise, a0_inb, sweep_half_c);

CL_alp_inb_M02 = tl.CL_alpha_DATCOM(AR_inb, 0.2, a0_inb, sweep_half_c);

CLmax_M02_inb = CLmax_Clmax*Clmax_M02_inb;				 # M = 0.2

alpha_stall_inb = np.degrees(CLmax_M02_inb/CL_alp_inb_M02 + np.radians(Alpha0_inb) + np.radians(d_alpha_clmax));

# OUTBOARD
CL_alp_outb = tl.CL_alpha_DATCOM(AR_outb, M_cruise, a0_outb, sweep_half_c);

CL_alp_outb_M02 = tl.CL_alpha_DATCOM(AR_outb, 0.2, a0_outb, sweep_half_c);

CLmax_M02_outb = CLmax_Clmax*Clmax_M02_outb;			 # M = 0.2

alpha_stall_outb = np.degrees(CLmax_M02_outb/CL_alp_outb_M02 + np.radians(Alpha0_outb) + np.radians(d_alpha_clmax));

# FULL WING
if alpha_stall_inb < alpha_stall_outb:
	alpha_stall = alpha_stall_inb;
elif alpha_stall_outb < alpha_stall_inb:
	alpha_stall = alpha_stall_outb;
else:
	print("Stall angle is no good ma boi");

CL_inb = np.radians(CL_alp_inb_M02)*alpha_stall;
CL_outb = np.radians(CL_alp_outb_M02)*alpha_stall;
#CL_max = (span_inboard*CL_inb + span_outboard*CL_outb)/span;
CL_max = (S_inb*CL_inb + S_outb*CL_outb)/S;


# ----- Drag Calculation -------
CD_misc = 5				# % of total CD0

# INBOARD
CD_0_cruise_inb = (1 + (CD_misc/100))*tl.CD0_wing(S_wet_ratio, t_c_avg_inb, x_c_m_inb, sweep_LE, C_r_m, C_t_m, 
					Cr_inb, span, rho_cruise, V_cruise, mu_cruise, k_wing, M_cruise, laminar_flow_inb, turb_flow_inb);

delta_CD_inb = tl.delta_wave_drag(M_cruise, t_c_stream_inb, sweep_LE, CL_cruise, M_crit_airfoil_inb);

CD_induced_inb = (CL_cruise**2)/(np.pi*AR_inb*e);

Total_CD_inb = CD_0_cruise_inb + delta_CD_inb + CD_induced_inb;

# OUTBOARD

CD_0_cruise_outb = (1 + (CD_misc/100))*tl.CD0_wing(S_wet_ratio, t_c_avg_outb, x_c_m_outb, sweep_LE, C_r_m, C_t_m, 
					Cr_outb, span, rho_cruise, V_cruise, mu_cruise, k_wing, M_cruise, laminar_flow_outb, turb_flow_outb);

delta_CD_outb = tl.delta_wave_drag(M_cruise, t_c_stream_outb, sweep_LE, CL_cruise, M_crit_airfoil_outb);

CD_induced_outb = (CL_cruise**2)/(np.pi*AR_outb*e);

Total_CD_outb = CD_0_cruise_outb + delta_CD_outb + CD_induced_outb;

t_final = time.time();

# ---- Print Results -----
print("Inboard Section Results ..... \n \n");
print("CL_alpha  = ", CL_alp_inb, "\n");
print("Total_CD  = ", Total_CD_inb  ,"\n");
print("CL_max    = ", CLmax_M02_inb,"  \n ");
print("Stall angle =", alpha_stall_inb,"   \n \n \n ");

print("Outboard Section Results ..... \n \n");
print("CL_alpha  =", CL_alp_outb,"\n");
print("Total_CD  =", Total_CD_outb,  "\n");
print("CL_max    =  ", CLmax_M02_outb," \n");
print("Stall angle = ", alpha_stall_outb,"  \n \n \n ");

print("Full Wing Results ..... \n \n");
print("Stall AOA of first section to stall =", alpha_stall, "\n");
print("Average CL_max at this AOA =", CL_max, "\n \n \n");

print("Optimal Taper Ratio:", planf.Optimal_taper(sweep_quart),"\n");

print("Time Taken:", t_final-t_s, "seconds");

print(sweep_half_c);
print(sweep_quart);


print(CL_alp_inb_M02);
print(CL_alp_outb_M02);

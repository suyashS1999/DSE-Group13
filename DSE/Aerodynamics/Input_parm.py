# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import Aerodynamics.WingPlanform as planf


#Global Wing

AR = 17;									# Wing Aspect ratio [-]
S = 150.8554095;							# Wing surface area [m^2]
span = 50.64130687;							# Wing span [m]
ave_chord = S/span;							# Average chord lenght [m]
taper_ratio = 0.44;							# Taper ratio [-]
# Inboard Wing
span_inboard = 36;							# Span of inboard wing [m]
span_outboard = span - span_inboard;		# Span of outboard wing [m]
t_c_avg  = 0.14         # Aver. t/c
x_c_m    = 0.1          # location of max thickness
# t_c_stream =           #-
CL_cruise  = 0.4370176576     #-

S_wet_ratio = 2.14         #- 

C_r_m, C_t_m = planf.Calc_root_tip_chordMain(span, AR, taper_ratio)

sweep_LE = np.radians(30) 

MAC = planf.Compute_MAC(C_r_m,C_t_m, sweep_LE, span)

#e  = e(AR,sweep_LE)
e = 0.8

rho_cruise = 0.3636     # kg/m^3
V_cruise   = 230.13     # m/s
mu_cruise = 1.46e-5     #SI units

laminar_flow = 0.1      #%
turb_flow = 1-laminar_flow  #%

k_wing = 0.634e-5        # Assuming smooth paint on wing

M_cruise = 0.78 


# Inboard

Inb_params, Outb_params = planf.InboardOutboard_wing_parms(S, span_inboard, span_outboard, C_r_m, C_t_m)

Cr_inb = Inb_params[0]
Ct_inb = Inb_params[1]
S_inb  = Inb_params[2]
AR_inb = Inb_params[3]
taper_ratio_inb = Inb_params[4]

a0_inb = 6.445775195   #1/rad
Cl0_inb = 0.314
Cl_max_inb = 2.26
Cd_inb  =0.006
alp_max_airf_inb = 22.5
t_c_avg_inb = 0.14
x_c_m_inb = 0.3
t_c_stream_inb = 0.14
M_crit_airfoil_inb = 0.5956992   #- 
laminar_flow_inb  = 0.1
turb_flow_inb = 0.9 

# Outboard

Cr_outb = Outb_params[0]
Ct_outb = Outb_params[1]
S_outb  = Outb_params[2]
AR_outb = Outb_params[3]
taper_ratio_outb = Outb_params[4]

a0_outb = 6.6864175   #1/rad
Cl0_outb = 0.45
Cl_max_outb = 1.6
Cd_outb  = 0.0045
alp_max_airf_outb = 20
t_c_avg_outb = 0.18
x_c_m_outb = 0.5
t_c_stream_outb = 0.18
M_crit_airfoil_outb = 0.62  #- 
laminar_flow_outb  = 0.1
turb_flow_outb = 0.9  


# If full wing Airfoil


a0_full = 0.10853430292060387*180/np.pi  # lift slope for airfoil (in 1/rad)
t_c_avg_full = 0.14 
x_c_m_full = 0.3636  # position of max thickness
laminar_flow_full = 0.1
turb_flow_full = 0.9

t_c_stream_full = 0.14
M_crit_airfoil_full = 0.5956992



# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import WingPlanform as planf


#Global Wing

AR = 17;									# Wing Aspect ratio [-]
S = 142.520611;								# Wing surface area [m^2]
span = 49.22245815;							# Wing span [m]
ave_chord = S/span;							# Average chord lenght [m]
taper_ratio = 0.44;							# Taper ratio [-]
# Inboard Wing
span_inboard = 36;							# Span of inboard wing [m]
span_outboard = span - span_inboard;		# Span of outboard wing [m]
t_c_avg  = 0.15         # Aver. t/c
x_c_m    = 0.1          # location of max thickness
# t_c_stream =           #-
CL_cruise  = 0.4744     #-

S_wet_ratio = 2.14         #- 

C_r_m, C_t_m = planf.Calc_root_tip_chordMain(span, AR, taper_ratio)

sweep_LE = np.radians(40.63) 

MAC = planf.Compute_MAC(C_r_m,C_t_m, sweep_LE, span)

# e  = e(AR,sweep_LE)
e = 0.7

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
Cl0_inb = 0.4
Cl_max_inb = 1.67
Cd_inb  = 0.005
alp_max_airf_inb = 15
t_c_avg_inb = 0.15
x_c_m_inb = 0.3
t_c_stream_inb = 0.15
M_crit_airfoil_inb = 0.592   #- 
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



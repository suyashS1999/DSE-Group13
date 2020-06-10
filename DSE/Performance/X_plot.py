#%% ---------------------- About this code ----------------------
# This code draws the scissor plot, made by the controllability curve and the stability curve
# It draws for both cruise situation and the approach situation


#%% ---------------------- Imports ----------------------
import matplotlib.pyplot as plt
import numpy as np

#%% ---------------------- Constants ----------------------


#%% ---------------------- Inputs ----------------------
#Geometric things
from Loading_diagram import cg_frw, cg_aft, xlemac
lmac = 1.501            # Mean aerodynamic chord [m]
#L0_angle = 0.1          # Zero lift angle, see the input Perf-23 for clarity [rad]
h_h = 3.6               # Height of the horizontal tail from the top of fuselage [m]
c_r = 4.021             # Root chord [m]
l_t = 42.6-6*0.75-1-xlemac+10.58-10.3                # Moment arm of horizontal tail [m]
b = 139.1               # Wing span [m]
Lamda_qc = 0.696        # Quarter chord sweep angle [rad]
#CL_alpha_w = 4.024      # (dummy) Cl_alpha of wing AND struts combined
AR = 17                 # Aspect ratio

#%% ---------------------- Functions ----------------------
def stablimit(ShS): 
    """
    Draws the stability curve, including the 5% safety margin
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app + (CL_alpha_h_app/CL_alpha_Ah_app)*(1 - downwash)*ShS*(l_t/lmac) - 0.05)
    xcg_cru = (xac_cru   + (CL_alpha_h_cru/CL_alpha_Ah_cru)*(1 - downwash)*ShS*(l_t/lmac) - 0.05)
        
    plt.plot(xcg_app, ShS, label = "Stability Approach")
    plt.plot(xcg_cru, ShS, label = "Stability Cruise")
    return (xcg_app, xcg_cru)

def contlimit(ShS):
    """
    Draws the controllability curve
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app - (Cm_ac_app/CL_Ah_app) + (CL_h/CL_Ah_app)*ShS*(l_t/lmac)*0.85)
    xcg_cru = (xac_cru - (Cm_ac_cru/CL_Ah_cru) + (CL_h/CL_Ah_cru)*ShS*(l_t/lmac)*0.85)
        
    plt.plot(xcg_app, ShS, label = "Controllability Approach")
    plt.plot(xcg_cru, ShS, label = "Controllability Cruise")
    return (xcg_app, xcg_cru)

#Functions below compute the parameters used by above plotting functions
#def downwash(zero_lift_angle):
#    mtv = np.cos(zero_lift_angle)*(h_h + (l_t-3/4*c_r)*np.tan(zero_lift_angle))*2/b
#    r = 2*l_t/b
#    K1 = (0.1124+0.1265*Lamda_qc+0.1766*Lamda_qc**2)/r**2+0.1024/r+2
#    K2 = 0.1124/r**2+0.1024/r+2
#    downwash = K1/K2*(r/(r**2+mtv**2)*0.4876/np.sqrt(r**2+0.6319+mtv**2)+(1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*CL_alpha_w/(np.pi*AR)
#    return downwash

xac_app_LERC = 10.263099  #position of ac during approach from the root chord leading edge [m]
xac_cru_LERC = 10.098687

ShS = np.linspace(0, 0.4, 2)
xac_app = (-10.58+xac_app_LERC)/lmac
xac_cru = (-10.58+xac_cru_LERC)/lmac
CL_alpha_h_app = 2.94719
CL_alpha_h_cru = 3.28587154
downwash = 0
CL_alpha_Ah_app = 0.065757*180/np.pi
CL_alpha_Ah_cru = 0.083393*180/np.pi

Cm_ac_app = -0.564400
Cm_ac_cru = -0.722830

CL_h = -0.05482782361203875

CL_Ah_app = 2 - CL_h
CL_Ah_cru = 0.566 - CL_h



stablimit(ShS)
contlimit(ShS)
frw = (cg_frw-xlemac)/lmac
aft = (cg_aft-xlemac)/lmac
plt.vlines(x = frw, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.vlines(x = aft, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.hlines(y = 0.26, xmin = frw, xmax = aft)
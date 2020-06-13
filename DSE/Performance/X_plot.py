#%% ---------------------- About this code ----------------------
# This code draws the scissor plot, made by the controllability curve and the stability curve
# It draws for both cruise situation and the approach situation


#%% ---------------------- Imports ----------------------
import matplotlib.pyplot as plt
import numpy as np

#%% ---------------------- Constants ----------------------
from Loading_diagram import cg_frw, cg_aft, xlemac, lmac
cg_aft = cg_aft
cg_frw = cg_frw
ymac = 11.01426343              #y position of mac from fuselage centre line [m]
Lamda_LE = 0.5235987756         #sweep angle LE [rad]
CL_cru = 0.5101828488               #CL during cruise [-]
CL_app = 2.4               #CL during approach [-]
S = 150.8519904                       #Surface area wing

#l_t = 42.6-6*0.75-1-xlemac+0.3 
AR = 17
h_h = 2.6
c_r = 4.135518182
b = 50.62
Lamda_qc = 0.4351399017
#%% ---------------------- Inputs ----------------------
#Aerodynamic parameters (from aero department)
xac_app_LERC = 6.446274                 # position of ac during approach from the root chord leading edge [m]
xac_cru_LERC = 6.859900                # position of ac during approach from the root chord leading edge [m]
CL_alpha_Ah_app = 0.066582*180/np.pi    # lift coeff gradient of aircraft-lesstail during approach [/rad]
CL_alpha_Ah_cru = 0.089685*180/np.pi    # lift coeff gradient of aircraft-lesstail during cruise [/rad]
CL_alpha_h_app = 2.94719                # lift coeff gradient of horizontal tail during approach [/rad]
CL_alpha_h_cru = 3.28587154             # lift coeff gradient of horizontal tail during approach [/rad]
Vh_V = 1                                # efficiency factor of H tail (or 1 if we opt for T-tail) [-]

Sh_S = 0.11
mu1 = 0.213
mu2 = 0.38
mu3 = 0.048
dCl = 1.9544
c_dash_c = 1.22
CL_TO = 2.1
Swf_S = 0.5

xlerc = xlemac - np.tan(Lamda_LE)*ymac     # X position of root chord LE
xac_app =  0.25 #(xlerc + xac_app_LERC-xlemac)/lmac - 0.398
xac_cru = (xlerc + xac_cru_LERC-xlemac)/lmac - 0.296
l_t = 42.6-6*0.75-1-(xlerc + xac_app_LERC)+3                # Moment arm of horizontal tail [m]

#%% ---------------------- Functions ----------------------
def stablimit(ShS): 
    """
    Draws the stability curve, including the 5% safety margin
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app + (CL_alpha_h_app/CL_alpha_Ah_app)*(1 - downwash)*ShS*(l_t/lmac)*Vh_V**2 - 0.05)
    xcg_cru = (xac_cru + (CL_alpha_h_cru/CL_alpha_Ah_cru)*(1 - downwash)*ShS*(l_t/lmac)*Vh_V**2 - 0.05)
        
    plt.plot(xcg_app, ShS, label = "Stability Approach")
    plt.plot(xcg_cru, ShS, label = "Stability Cruise")
    return xcg_app, xcg_cru

def contlimit(ShS):
    """
    Draws the controllability curve
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app - (Cm_ac_app/CL_Ah_app) + (CL_h_app/CL_Ah_app)*ShS*(l_t/lmac)*Vh_V**2)
    xcg_cru = (xac_cru - (Cm_ac_cru/CL_Ah_cru) + (CL_h_cru/CL_Ah_cru)*ShS*(l_t/lmac)*Vh_V**2)
        
    plt.plot(xcg_app, ShS, label = "Controllability Approach")
#    plt.plot(xcg_cru, ShS, label = "Controllability Cruise")
    return xcg_app, xcg_cru

def req_CLh():
    """
    determines required CL_h for a given CG range and the CL required
    """
    Cm_ac = [Cm_ac_app, Cm_ac_cru]
    CL = [CL_app, CL_cru]
    CL_h_app = (-Cm_ac[0] + CL[0]*(cg_frw-(xlerc + xac_app_LERC)))/(Vh_V**2*Sh_S*l_t/lmac), (-Cm_ac[0] + CL[0]*(cg_aft-(xlerc + xac_app_LERC)))/(Sh_S*l_t/lmac)
    CL_h_cru = (-Cm_ac[1] + CL[1]*(cg_frw-(xlerc + xac_cru_LERC)))/(Vh_V**2*Sh_S*l_t/lmac), (-Cm_ac[1] + CL[1]*(cg_aft-(xlerc + xac_cru_LERC)))/(Sh_S*l_t/lmac)
    return CL_h_app, CL_h_cru

def downwash_tail(zero_lift_angle):
    mtv = np.cos(zero_lift_angle)*(h_h + (l_t-3/4*c_r)*np.tan(zero_lift_angle))*2/b
    r = 2*l_t/b
    K1 = (0.1124+0.1265*Lamda_qc+0.1766*Lamda_qc**2)/r**2+0.1024/r+2
    K2 = 0.1124/r**2+0.1024/r+2
    downwash = K1/K2*(r/(r**2+mtv**2)*0.4876/np.sqrt(r**2+0.6319+mtv**2)+(1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*CL_alpha_Ah_cru/(np.pi*AR)
    return downwash

def Cm_flaps():
    dCM_quart = mu2*(-mu1*dCl*c_dash_c - (CL_TO + dCl*(1-Swf_S))*1/8*c_dash_c*(c_dash_c - 1)) + 0.7 * AR/(1+2/AR) *mu3 *dCl*np.tan(Lamda_qc)
    CM_ac_flap = dCM_quart - CL_app * (0.25 - xac_app)
    return CM_ac_flap

#%% ---------------------- Main ----------------------
Cm_ac_app = +0.177053 - Cm_flaps()                   # pitching moment coefficient at aero-centre during approach [-]
Cm_ac_cru = +0.278472                   # pitching moment coefficient at aero-centre during approach [-]
#counter clockwise positive!

#downwash = downwash_tail(2)                            #No disturbance from wing downwash (d epsilon/ d alpha = 0) [-]
downwash = 0

CL_h_app, CL_h_cru = -0.505, req_CLh()[1][1]
#CL_h_app, CL_h_cru = max(CL_h_app, key=abs), max(CL_h_cru, key=abs)


CL_Ah_app = CL_app - CL_h_app
CL_Ah_cru = CL_cru - CL_h_cru

Cm_a = CL_alpha_Ah_app * (frw - xac_app) - CL_alpha_h_app*Sh_S*l_t/lmac

#Making the plot below
ShS = np.linspace(0, 0.4, 2)
stablimit(ShS)
contlimit(ShS)
plt.legend()
frw = (cg_frw-xlemac)/lmac
aft = (cg_aft-xlemac)/lmac
plt.vlines(x = frw, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.vlines(x = aft, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.hlines(y = Sh_S, xmin = frw, xmax = aft)
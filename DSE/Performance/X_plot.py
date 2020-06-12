#%% ---------------------- About this code ----------------------
# This code draws the scissor plot, made by the controllability curve and the stability curve
# It draws for both cruise situation and the approach situation


#%% ---------------------- Imports ----------------------
import matplotlib.pyplot as plt
import numpy as np

#%% ---------------------- Constants ----------------------
from Loading_diagram import cg_frw, cg_aft, xlemac, lmac
ymac = 11.01426343              #y position of mac from fuselage centre line [m]
Lamda_LE = 0.5235987756         #sweep angle LE [rad]
CL_cru = 0.566               #CL during cruise [-]
CL_app = 2.0               #CL during approach [-]
S = 150.7                       #Surface area wing
#%% ---------------------- Inputs ----------------------
#Aerodynamic parameters (from aero department)
xac_app_LERC = 10.263099                #position of ac during approach from the root chord leading edge [m]
xac_cru_LERC = 10.098687                #position of ac during approach from the root chord leading edge [m]
CL_alpha_Ah_app = 0.065757*180/np.pi    #lift coeff gradient of aircraft-lesstail during approach [/rad]
CL_alpha_Ah_cru = 0.083393*180/np.pi    #lift coeff gradient of aircraft-lesstail during cruise [/rad]
Cm_ac_app = -0.564400                   #pitching moment coefficient at aero-centre during approach [-]
Cm_ac_cru = -0.722830                   #pitching moment coefficient at aero-centre during approach [-]
CL_alpha_h_app = 2.94719                #lift coeff gradient of horizontal tail during approach [/rad]
CL_alpha_h_cru = 3.28587154             #lift coeff gradient of horizontal tail during approach [/rad]
downwash = 0                            #No disturbance from wing downwash (d epsilon/ d alpha = 0) [-]
Vh_V = 0.95                             #efficiency factor of H tail (or 1 if we opt for T-tail) [-]

Sh_S = 0.25

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
    return (xcg_app, xcg_cru)

def contlimit(ShS):
    """
    Draws the controllability curve
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app - (Cm_ac_app/CL_Ah_app) + (CL_h/CL_Ah_app)*ShS*(l_t/lmac))
    xcg_cru = (xac_cru - (Cm_ac_cru/CL_Ah_cru) + (CL_h/CL_Ah_cru)*ShS*(l_t/lmac))
        
    plt.plot(xcg_app, ShS, label = "Controllability Approach")
    plt.plot(xcg_cru, ShS, label = "Controllability Cruise")
    return (xcg_app, xcg_cru)

def req_CLh():
    """
    determines required CL_h for a given CG range and the CL required
    """
    Cm_ac = [Cm_ac_app, Cm_ac_cru]
    CL = [CL_app, CL_cru]
    CL_h_app = (Cm_ac[0] + CL[0]*(cg_frw-xac_app))/(Sh_S*l_t/lmac), (Cm_ac[0] + CL[0]*(cg_aft-xac_app))/(Sh_S*l_t/lmac)
    CL_h_cru = (Cm_ac[1] + CL[1]*(cg_frw-xac_app))/(Sh_S*l_t/lmac), (Cm_ac[1] + CL[1]*(cg_aft-xac_app))/(Sh_S*l_t/lmac)
    return CL_h_app, CL_h_cru

#%% ---------------------- Main ----------------------


xlerc = xlemac - np.tan(Lamda_LE)*ymac     # X position of root chord LE
xac_app = (xlerc + xac_app_LERC-xlemac)/lmac
xac_cru = (xlerc + xac_cru_LERC-xlemac)/lmac
l_t = 42.6-6*0.75-1.2-xac_app                # Moment arm of horizontal tail [m]

CL_h_app, CL_h_cru = req_CLh()

CL_Ah_app = CL_app - CL_h_app[1]        #toggle between [0] and [1] to see which is more critical
CL_Ah_cru = CL_cru - CL_h_cru[1]

#Making the plot below
ShS = np.linspace(0, 0.4, 2)
stablimit(ShS)
contlimit(ShS)
frw = (cg_frw-xlemac)/lmac
aft = (cg_aft-xlemac)/lmac
plt.vlines(x = frw, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.vlines(x = aft, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.hlines(y = Sh_S, xmin = frw, xmax = aft)
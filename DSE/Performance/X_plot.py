#%% ---------------------- About this code ----------------------
# This code draws the scissor plot, made by the controllability curve and the stability curve
# It draws for both cruise situation and the approach situation


#%% ---------------------- Imports ----------------------
import matplotlib.pyplot as plt
import numpy as np

#%% ---------------------- Constants ----------------------
from Loading_diagram import cg_frw, cg_aft, xlemac, lmac, cg_frw_mac, cg_aft_mac
frw = cg_frw_mac
aft = cg_aft_mac
ymac = 11.01426343              #y position of mac from fuselage centre line [m]
Lamda_LE = 0.5235987756         #sweep angle LE [rad]
CL_cru = 0.4890523925               #CL during cruise [-]
CL_app = 2.4               #CL during approach [-]
S = 150.7                       #Surface area wing

#l_t = 42.6-6*0.75-1-xlemac+0.3 
AR = 17
h_h = 1.2
c_r = 4.135518182
b = 50.62
Lamda_qc = 0.4351399017     #sweep angle quarter chord [rad]
#%% ---------------------- Inputs ----------------------
#Aerodynamic parameters (from aero department)
#xac_app_LERC = 6.859900                  # position of ac during approach from the root chord leading edge [m] VSP #Suyash check please
#xac_cru_LERC = 6.446274     #now its switched           # position of ac during approach from the root chord leading edge [m] VSP
CL_alpha_Ah_app = 0.066582*180/np.pi    # lift coeff gradient of aircraft-lesstail during approach [/rad] VSP
CL_alpha_Ah_cru = 0.089685*180/np.pi    # lift coeff gradient of aircraft-lesstail during cruise [/rad] VSP
CL_alpha_h_app = 2.94719                # lift coeff gradient of horizontal tail during approach [/rad] DATCOM
CL_alpha_h_cru = 3.28587154             # lift coeff gradient of horizontal tail during approach [/rad] DATCOM
Vh_V = np.sqrt(0.95)                             # efficiency factor of H tail (or 1 if we opt for T-tail) [-]

Sh_S = 0.12
mu1 = 0.213
mu2 = 0.38
mu3 = 0.048
dCl = 1.9544 #?
c_dash_c = 1.22     #flap extention ratio
CL_TO = 2.1
Swf_S = 0.5

xlerc = xlemac - np.tan(Lamda_LE)*ymac     # X position of root chord LE
l_fn = 1.212 + xlerc
b_f = 3.9
h_f = 4.27
lamda = 0.44
x_ac_w_app = 0.3
x_ac_w_cru = 0.4
x_ac_f_app = (-1.8/CL_alpha_Ah_app)*b_f*h_f*l_fn/(S*lmac) + 0.273/(1+lamda)*b_f*S/b*(b-b_f)/(lmac**2*(b+2.15*b_f))*np.tan(Lamda_qc)
x_ac_f_cru = (-1.8/CL_alpha_Ah_cru)*b_f*h_f*l_fn/(S*lmac) + 0.273/(1+lamda)*b_f*S/b*(b-b_f)/(lmac**2*(b+2.15*b_f))*np.tan(Lamda_qc)
xac_app = x_ac_w_app + x_ac_f_app - 0.554
xac_cru = x_ac_w_cru + x_ac_f_cru - 0.411

l_t_app = 42.6-6*0.75-1-(xlemac + xac_app*lmac)+1                # Moment arm of horizontal tail approach [m]
l_t_cru = 42.6-6*0.75-1-(xlemac + xac_cru*lmac)+1   # Moment arm of horizontal tail cruise [m]

#xac_app = (xlerc + xac_app_LERC-xlemac)/lmac - 0.554 #Converting to MAC
#xac_cru = (xlerc + xac_cru_LERC-xlemac)/lmac - 0.411 #Converting to MAC

#%% ---------------------- Functions ----------------------
def stablimit(ShS): 
    """
    Draws the stability curve, including the 5% safety margin
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app + (CL_alpha_h_app/CL_alpha_Ah_app)*(1 - downwash)*ShS*(l_t_app/lmac)*Vh_V**2 - 0.05)
    xcg_cru = (xac_cru + (CL_alpha_h_cru/CL_alpha_Ah_cru)*(1 - downwash)*ShS*(l_t_cru/lmac)*Vh_V**2 - 0.05)
        
    plt.plot(xcg_app, ShS, label = "Stability Approach")
    plt.plot(xcg_cru, ShS, label = "Stability Cruise")
    return xcg_app, xcg_cru

def contlimit(ShS):
    """
    Draws the controllability curve
    Only input is the ShS, which should be a numpy array with length of 2, like [0, 0.5]
    Other inputs are implicitly taken, which should all be computed before this function is ran
    """
    xcg_app = (xac_app - (Cm_ac_app/CL_Ah_app) + (CL_h_app/CL_Ah_app)*ShS*(l_t_app/lmac)*Vh_V**2)
    xcg_cru = (xac_cru - (Cm_ac_cru/CL_Ah_cru) + (CL_h_cru/CL_Ah_cru)*ShS*(l_t_cru/lmac)*Vh_V**2)
        
    plt.plot(xcg_app, ShS, label = "Controllability Approach")
    plt.plot(xcg_cru, ShS, label = "Controllability Cruise")
    return xcg_app, xcg_cru

def req_CLh():
    """
    determines required CL_h for a given CG range and the CL required
    """
    Cm_ac = [Cm_ac_app, Cm_ac_cru]
    CL = [CL_app, CL_cru]
    CL_h_app = (Cm_ac[0] + CL[0]*(frw-(xac_app)))/(Vh_V**2*Sh_S*l_t_app/lmac), (Cm_ac[0] + CL[0]*(aft-xac_app))/(Sh_S*l_t_app/lmac)
    CL_h_cru = (Cm_ac[1] + CL[1]*(frw-(xac_cru)))/(Vh_V**2*Sh_S*l_t_cru/lmac), (Cm_ac[1] + CL[1]*(aft-xac_cru))/(Sh_S*l_t_cru/lmac)
    return CL_h_app, CL_h_cru

def downwash_tail(zero_lift_angle):
    mtv = np.cos(zero_lift_angle)*(h_h + (l_t_app-3/4*c_r)*np.tan(zero_lift_angle))*2/b
    r = 2*l_t_app/b
    K1 = (0.1124+0.1265*Lamda_qc+0.1766*Lamda_qc**2)/r**2+0.1024/r+2
    K2 = 0.1124/r**2+0.1024/r+2
    downwash = K1/K2*(r/(r**2+mtv**2)*0.4876/np.sqrt(r**2+0.6319+mtv**2)+(1+(r**2/(r**2+0.7915+5.0734*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*CL_alpha_Ah_cru/(np.pi*AR)
    return downwash

def Cm_flaps():
    dCM_quart = mu2*(-mu1*dCl*c_dash_c - (CL_TO + dCl*(1-Swf_S))*1/8*c_dash_c*(c_dash_c - 1)) + 0.7 * AR/(1+2/AR) *mu3 *dCl*np.tan(Lamda_qc)
    CM_ac_flap = dCM_quart - CL_app * (0.25 - xac_app)
    CM_ac_flap = -0.19
    return CM_ac_flap

#%% ---------------------- Main ----------------------
#Cm_ac_app = +0.177053 - Cm_flaps()                   # pitching moment coefficient at aero-centre during approach [-]
#Cm_ac_cru = -0.278472                   # pitching moment coefficient at aero-centre during approach [-]
#counter clockwise positive!
Cm_ac_app = -0.0746 - 0.2457 + Cm_flaps()   #wing contribution + fuselage contribution + flaps contribution to Cm around ac during cruise
Cm_ac_cru = -0.0746 - 0.1148                #wing contribution + fuselage contribution to Cm around ac during cruise


downwash = downwash_tail(2.7)                            #No disturbance from wing downwash (d epsilon/ d alpha = 0) [-]
#downwash = 0

CL_h_app = -0.505 #maximum negative lift h tail can deliver (dependant on AR, from SEAD)
CL_h_cru = -0.505
#CL_h_app, CL_h_cru = max(CL_h_app, key=abs), max(CL_h_cru, key=abs)

l_fn = 1.212 + xlerc
b_f = 3.92
h_f = 4.37
lamda = 0.44
x_ac_w_app = 0.35
x_ac_w_cru = 0.4
x_ac_f_app = (-1.8/CL_alpha_Ah_app)*b_f*h_f*l_fn/(S*lmac) + 0.273/(1+lamda)*b_f*S/b*(b-b_f)/(lmac**2*(b+2.15*b_f))*np.tan(Lamda_qc)
x_ac_f_cru = (-1.8/CL_alpha_Ah_cru)*b_f*h_f*l_fn/(S*lmac) + 0.273/(1+lamda)*b_f*S/b*(b-b_f)/(lmac**2*(b+2.15*b_f))*np.tan(Lamda_qc)
xac_app = x_ac_w_app + x_ac_f_app - 0.554
xac_cru = x_ac_w_cru + x_ac_f_cru - 0.411

l_t_app = 42.6-6*0.75-1-(xlemac + xac_app*lmac)+1                # Moment arm of horizontal tail approach [m]
l_t_cru = 42.6-6*0.75-1-(xlemac + xac_cru*lmac)+1   # Moment arm of horizontal tail cruise [m]

CL_Ah_app = CL_app-CL_h_app
CL_Ah_cru = CL_cru*1.1
CL_h_cru_ave = np.average(req_CLh()[1])
#Cm_a = CL_alpha_Ah_app * (frw - xac_app) - CL_alpha_h_app*Sh_S*l_t/lmac

#Making the plot below
ShS = np.linspace(0, 0.4, 2)
stablimit(ShS)
contlimit(ShS)
plt.legend(fontsize = 12, loc = 'upper left')
plt.grid()
plt.xlabel("CG location in X direction [MAC]", fontsize = 16)
plt.ylabel("Sh/S [-]", fontsize = 16)
plt.vlines(x = frw, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.vlines(x = aft, ymin = 0, ymax = ShS[-1], linestyles = "--")
plt.hlines(y = Sh_S, xmin = frw, xmax = aft)
plt.xlim(-1, 2.5)
plt.ylim(0, 0.6)
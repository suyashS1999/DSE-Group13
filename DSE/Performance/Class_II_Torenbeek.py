#About this code
"""
This code performs the class II weight estimation with the method by Torenbeek
Everything in this code will be computed in SI units
"""

import numpy as np

#%% ---------------------- Constants ----------------------
k_s = 0.447     #Correction factor, airframe structure
k_w = 6.67e-3   #Correction factor, wing
b_ref = 1.905   #Span constant for wing sizing
k_h = 1.0       #Correction factor, h tail (fixed stabiliser)
k_wf = 0.23     #Correction factor, fuselage
k_uc = 1.08     #Correction factor, undercarriage due to high wing configuration
k_sc = 0.64*(1+0.2+0.15)*0.768     #Correction factor, undercarriage (with LE flap or slat controls + lift dumper controls assumed) (0.768 to correct for SI units)
k_pg = 1.15     #Correction factor, propulsion group, as podded-engine jet transports
k_thr = 1.0     #Correction factor, thrust (no thrust reverser hence 1.00)

#%% ------------------------ INPUTS ------------------------
#Inputs
n_ult = 3.75        #ultimate load factor [-]
MTOW = 74616.9829   #[kg]
MF = 12471.21232    #Fuel mass [kg]

V_D = 156.7297      #dive speed in EAS [m/s]
l_t = 12              #Length between quarter chords of wing and h tail [m]

b_f = 4.2              #width fuselage [m]
h_f = 4.2             #height fuselage [m]
l_f = 40            #length fuselage [m]
S_G = 450           #gross shell area (entire outer surface of the fuselage) [m^2]

b =  49.22               #wing span [m]
Lamda_halfc =   0.682    #half chord sweep [rad]
S =     142.5            #wing area [m^2]
t_r =    0.18*4.021           #maximum thickness root chord [m]

S_h =    30           #horizontal tail area [m^s]
Lamda_h = 0.3          #Sweep angle of horizontal tail [rad]
h_h =    6           #height of horizontal stabiliser
S_v =    20           #vertical tail area [m^s]
Lamda_v =  0.3         #Sweep angle of vertical tail [rad] 
b_v = h_h           #span of vertical tail (assumed to be same as height of the h tail because its a T tail)

T_TO = 211545.862   #Take off thrust [N] (From wing/thrust loading)
N_e = 2             #number of engine [-]
W_engine =  2500        #Total engine weight [kg]

#%% ------------------------ Functions ------------------------ 
#Functions weight estimations
#Airframe structure
def W_struct(): #Torenbeek eq. (8-11)
    W_s = MTOW*k_s*np.sqrt(n_ult)*(b_f*h_f*l_f/MTOW)**0.24
    return W_s
#Wing group 
def W_wing(): #Torenbeek eq. (8-12) (wrong, MTOW should be MZFW)
    MZFW = MTOW - MF    #Maximum zero fuel weight [kg]
    b_s = b/np.cos(Lamda_halfc)
    W_w = MZFW*k_w*b_s**0.75*(1+np.sqrt(b_ref/b_s))*n_ult**0.55*((b_s/t_r)/(MZFW/S))**0.3
    W_w = W_w*0.7 #Accounting for braced wing, including the struts
    return W_w #Verified

#Tail group
def W_tail():
    P_h = S_h**0.2*V_D/1000/np.sqrt(np.cos(Lamda_h))
    if P_h < 0.45: #If in linear part, compute automatically
        Q_h = 61.73*P_h-3.02
    elif P_h > 0.65:
        Q_h = 29.47
    else:
        print("The SV/cos is outside the linear section, it is", P_h, "(in the top axis with SI units)")
        Q_h = float(input("Look at the Torenbeek tail graph in this folder, input the y value (Use right axis)"))
    W_h = S_h*k_h*Q_h #Q_h and Q_v read off the graph
    P_v = S_v**0.2*V_D/1000/np.sqrt(np.cos(Lamda_v))
    if P_v < 0.45: #If in linear part, compute automatically
        Q_v = 61.73*P_h-3.02
    elif P_v > 0.65:
        Q_v = 29.47
    else:
        print("The SV/cos is outside the linear section, it is", P_v, "(in the top axis with SI units)")
        Q_v = float(input("Look at the Torenbeek tail graph in this folder, input the y value (Use right axis)"))
    k_v = 1+0.15*(S_h*h_h)/(S_v*b_v)    #Correction factor, for T-tail
    W_v = S_v*k_v*Q_v #Q_h and Q_v read off the graph
    return W_h + W_v

#Body group
def W_body(): #Torenbeek eq (8-16)
    W_f = k_wf * np.sqrt(V_D*l_t/(b_f+h_f))*S_G**1.2
    W_f = W_f+0.08*W_f+0.07*W_f #Weight increase from pressurised cabin and fuselage mounted landing gear, respectively
    return W_f #Verified

#Undercarriage
def W_undercarriage(): #eq (8-17)
    A_main = 18.1
    B_main = 0.131
    C_main = 0.019
    D_main = 2.23e-5
    W_uc_main = k_uc*(A_main+B_main*MTOW**(3/4)+C_main*MTOW+D_main*MTOW**(3/2))
    A_nose = 9.1
    B_nose = 0.082
    C_nose = 0
    D_nose = 2.97e-6
    W_uc_nose = k_uc*(A_nose+B_nose*MTOW**(3/4)+C_nose*MTOW+D_nose*MTOW**(3/2))
    return W_uc_main+W_uc_nose #Verified

#Surface controls group
def W_surfcont(): # eq (8-18)
    W_sc = k_sc*MTOW**(2/3)
    return W_sc #Verified

#Engine section and nacelle group
def W_nacelle(): # eq (8-25)
    W_n = 0.065*T_TO/9.81*0.9
    return W_n #Verified

#Propulsion group
def W_prop(): #eq (8-27)
    W_pg = k_pg*k_thr*N_e*W_engine
    return W_pg #Verified

#Airframe services and equipment
def W_equip():
    W_se= 0.11*MTOW #assume "medium range transport aircraft"
    return W_se #Verified

#%% ------------------------ Main ------------------------
W_airframe = W_wing() + W_tail() + W_body() + W_undercarriage() + W_surfcont() + W_nacelle()
W_prop = W_prop()
W_equip = W_equip()
OEW = W_airframe + W_prop + W_equip
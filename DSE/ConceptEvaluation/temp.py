# -*- coding: utf-8 -*-
"""
Used in the conceptual stage to estimate main parameters for concept selection
Outputs are
Payload range diagram
Weight
Volume
Cost 
"""

from numpy import*



# ------------------------ INPUTS ------------------------

# ---------------------- Constants ----------------------
g = 9.81
M = 0.78


# ---------------------- Variables ----------------------
h = 10970 #km = FL360

cj = 1.132*10**(-5)
L_D = 13
Weight_frac = 1/0.8544

# ------------------------ RANGE ------------------------



def ISA_trop(h):
    """ This function computes
	Input:
		
	Output:
		
	"""
    T = 288.15 -0.0065*h
    p = 101325*(T/288.15)**(-g/(-0.0065*287))
    rho = 1.225 * (T/288.15)**(-g/(-0.0065*287)-1)
    a = sqrt(1.4*287*T)
    
    return T,p,rho,a
    



def RangeJet(g, M, a, cj, L_D, Weight_frac):
    
    """ This function computes the range that can be obtained by a 
		jet aircraft
	Input:
		V = cruise velocity [m/s]
		cj = Thrust specific fuel consumption [kg/Ns]
		L_D = Lift to Drag ratio
		Weight_frac = weight fraction before and after cruise phase [-]
	Output:
		R = Range [m]
	"""
    
    V = M * a
    R = (V/(cj*g))*L_D*log(Weight_frac)
    
    return R



# ------------------------ MAIN -----------------------

T, p, rho, a = ISA_trop(h)
R = RangeJet(g, M, a, cj, L_D, Weight_frac)

print(R)

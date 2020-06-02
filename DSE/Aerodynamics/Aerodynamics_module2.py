import numpy as np



def Re_subsonic(rho,V,l,mu,k):
    R_actual = rho*V*l/mu
    R_cutoff = 38.21*((l/k)**1.053)
    
    return min(R_actual,R_cutoff)


def Re_transonic(rho,V,l,mu,k,M):
     R_actual = rho*V*l/mu
     
     R_cutoff = 44.62*((l/k)**1.053)*(M**1.16)
     
     return min(R_actual,R_cutoff)


def Cf_laminar(Re):
    return 1.328/np.sqrt(Re)

def Cf_turb(M,Re):
    
    cf = 0.455/(((np.log10(Re))**2.58)*(1+(0.144*(M**2)))**0.65)
    
    return cf


def Form_factor_wing(t_c,x_c_m,M,lmb_m):
    
    FF = (1+((0.6/x_c_m)*t_c)+(100*(t_c**4)))*(1.34*(M**0.18)*(np.cos(lmb_m)**0.28))
    
    return FF


def M_dd(t_c_stream,lmb_le,CL):
    
    Md = (0.935/np.cos(lmb_le)) - (t_c_stream/(np.cos(lmb_le)**2)) - (CL/(10*(np.cos(lmb_le)**3)))
    
    return Md


def sweep_x(x,lmb_le,Cr,Ct,b_total):
    
    sweep_arb = np.arctan((np.tan(lmb_le)) - (Cr*x/(b_total*0.5)) + (Ct*x/(b_total*0.5)))
    
    return sweep_arb
    

def e(AR,sweep_le):
    
    e = (4.61*(1-(0.045*(AR**0.68)))*((np.cos(sweep_le))**0.15)) - 3.1
    
    return e



import numpy as np

pe = 2e-5

def PBL(Pf,dfreq,pe0):

    PBL = 10*np.log10((Pf*dfreq)/(pe0**2))
    
    return PBL

def SPL(pe,pe0):

    SPL = 10*np.log10((pe**2)/(pe0**2))

    return SPL

def dL_A(f):

    dL_A = -145.528 + 98.262*np.log10(f) - 19.509*(np.log10(f))**2 + 0.975*(np.log10(f))**3

    return dL_A

def L_A(SPL,dL_A):
    
    L_A = 10*np.log10(np.sum(10**((SPL+dL_A)/10)))

    return L_A


import numpy as np

def compressor(T1,p1,eff,PI,mf,cp,K):
    """
    Compressor
    This functions outputs total temperature and pressure, and powered required for compressor.
    Inputs: total temp, total pres, efficiency comp, pressure ratio comp, massflow air, cp(air),K(air)
    """
    T2 = T1*(1+1/eff*(PI**((K-1)/K)-1))
    p2 = PI*p1
    Power = mf*cp*(T2-T1)
    return T2,p2,Power

def combustion(T1,T2,mf,LHV,cp,eff,p1,PI):
    """
    Combustion chamber
    This functions outputs mass flow fuel, total pressure, and massflow air + fuel.
    Inputs: total temp before, total temp after, massflow air + fuel, LHV_f, cp(gas),
            efficiency combustion, total pressure before, pressure ratio combustion
    """
    mf_f = mf*cp*eff*(T2-T1)/LHV
    p2 = p1*PI
    new_mf = mf + mf_f
    return mf_f, p2, new_mf

def turbine(Power,mf,cp,T1,p1,eff,K):
    """
    Turbine
    This functions outputs total temperature and pressure, and turbine pressure ratio.
    Inputs: power to deliver, massflow air + fuel, cp(gas), total temperature before
            total pressure before, efficiency turbine, K(air)
    """
    T2 = T1-Power/(mf*cp)
    p2 = p1*(1-1/eff*(1-T2/T1))**(K/(K-1))
    PI = p1/p2
    return T2, p2, PI

def exhaust(p0,T01,p01,eff,K,cp):
    """
    Exhaust
    This functions outputs temperature exhaust and velocity exhaust.
    Inputs: freestream pressure, total temperature before, total pressure before,
            efficiency of the exhaust, K(gas), cp(gas)
    """
    T_e = T01*(1-eff*(1-p0/p01))**((K-1)/K)
    V_e = np.sqrt(2*cp*(T01-T_e))
    return T_e, V_e

def thrust_core(T_e,T_0,mf):
    """
    Thrust from core
    This functions outputs core thrust.
    Inputs: temperature exit, temperature before, massflow core fuel + air
    """
    Thr_core = mf*(T_e-T_0) #assumed no pressure difference
    return Thr_core
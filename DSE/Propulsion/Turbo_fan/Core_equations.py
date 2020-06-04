import numpy as np

def combustion(T1,T2,mf,LHV,cp,eff,p1,PI):
    #Total temperatures and pressures used
    #use efficiency of the combustion chamber
    #use hot cp and hot K
    #massflow without fuel
    #PI is pressure ratio combustion chamber
    mf_f = mf*cp*eff*(T2-T1)/LHV
    p2 = p1*PI
    new_mf = mf + mf_f
    return mf_f, p2, new_mf
    #p2 is total temperature after combustion
    #mf_f is the massflow of the fuel
    #new_mf is the new massflow out of the combustion chamber

def turbine(Power,mf,cp,T1,p1,eff,K):
    #Total temperatures and pressures used
    #use efficiency of the turbine
    #use hot cp and hot K
    #Power is the power the turbine has to provide
    T2 = T1-Power/(mf*cp)
    p2 = p1*(1-1/eff*(1-T2/T1))**(K/(K-1))
    PI = p1/p2
    return T2, p2, PI
    #T2 is total temp after turbine
    #p2 is total pressure after turbine
    #PI is pressure ratio of the turbine

def exhaust(p0,T01,eff,p01,K,cp):
    #p0 is pressure in the freestream, not total pressure!!!
    #T01 is total temperature after turbine
    #eff is efficiency of the nozzle
    #p01 is total pressure after turbine
    #K and cp hot
    T_e = T01*(1-eff*(1-p0/p01))**((K-1)/K)
    V_e = np.sqrt(2*cp*(T01-T_e))
    return T_e, V_e
    #returns exit temperature and exit velocity


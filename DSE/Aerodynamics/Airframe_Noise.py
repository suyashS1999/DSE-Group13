import numpy as np
from Noise_tools import *
from HLDsizing import *
from Input_parm import *

pe0 = 2e-5

f = np.array([1.25,1.6,2,2.5,3.15,4,5,6.3,8,10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000])

#Inputs
Aw = 150.8519904                  #Area of the wing [m2]
bw = 50.64073298                  #Span of the wing [m]
rho_inf = 1.225             #Rho at cruise [m3/kg]
M = 0.2                     #Mach number [-]
c = 340.3                   #Speed of sound [m/s]
dynamic_visc = 1.84e-5      #Dynamic viscosity in the Netherlands at Schiphol[kg/(ms)]
flapped_area = 0.50 #percent
Af = 0.25*flapped_area*Aw                     #Total flap area [m2]
bf = 23                       #Total flap span [m]
df = HLD.TE_HLDs[TE_HLD][1] #Flap deflection angle [deg]
d_main = 1.1684                       #Wheel diameter [m]
d_nose = 0.762                       #Wheel diameter [m]
n_main = 2                  #Number of wheels [-]
n_nose = 2                  #Number of wheels [-]

phi = 0                         #Azimuthal directivity angle [deg]
r = 1                           #distance from source to observer
thetal = np.linspace(0,360,361)  #Polar directivity angle [deg]

#Start calculations
pe_sq = np.zeros((len(f),len(thetal)))

for i in range(len(thetal)):
    theta = thetal[i]
    theta = np.array([theta])

    G_cleanwing = geometry_func(Aw,bw,rho_inf,M,c,dynamic_visc,df,n_main,d_main,"Clean wing")
    G_slats = geometry_func(Aw,bw,rho_inf,M,c,dynamic_visc,df,n_main,d_main,"Slats")
    G_flaps = geometry_func(Af,bw,rho_inf,M,c,dynamic_visc,df,n_main,d_main,"Flaps")
    G_mainlg = geometry_func(Aw,bw,rho_inf,M,c,dynamic_visc,df,n_main,d_main,"Main landing gear")
    G_noselg = geometry_func(Aw,bw,rho_inf,M,c,dynamic_visc,df,n_nose,d_nose,"Main landing gear")

    L_cleanwing = length_scale(G_cleanwing,bw,Aw,d_main,"Clean wing")
    L_slats = length_scale(G_slats,bw,Aw,d_main,"Slats")
    L_flaps = length_scale(G_flaps,bf,Af,d_main,"Flaps")
    L_mainlg = length_scale(G_mainlg,bf,Af,d_main,"Main landing gear")
    L_noselg = length_scale(G_noselg,bf,Af,d_nose,"Main landing gear")

    K_cleanwing = constants_Ka("Clean wing",n_main)[0]
    a_cleanwing = constants_Ka("Clean wing",n_main)[1]
    K_slats = constants_Ka("Slats",n_main)[0]
    a_slats = constants_Ka("Slats",n_main)[1]
    K_flaps = constants_Ka("Flaps",n_main)[0]
    a_flaps = constants_Ka("Flaps",n_main)[1]
    K_mainlg = constants_Ka("Main landing gear",n_main)[0]
    a_mainlg = constants_Ka("Main landing gear",n_main)[1]
    K_noselg = constants_Ka("Main landing gear",n_nose)[0]
    a_noselg = constants_Ka("Main landing gear",n_nose)[1]


    P_cleanwing = powers(K_cleanwing, M, a_cleanwing, G_cleanwing,rho_inf,c,bw)
    P_slats = powers(K_slats,M,a_slats,G_slats,rho_inf,c,bw)
    P_flaps = powers(K_flaps,M,a_flaps,G_flaps,rho_inf,c,bw)
    P_mainlg = powers(K_mainlg,M,a_mainlg,G_mainlg,rho_inf,c,bw)
    P_noselg = powers(K_noselg,M,a_noselg,G_noselg,rho_inf,c,bw)

    S_cleanwing = Strouhal_number(f,L_cleanwing,M,theta,c)
    S_slats = Strouhal_number(f,L_slats,M,theta,c)
    S_flaps = Strouhal_number(f,L_flaps,M,theta,c)
    S_mainlg = Strouhal_number(f,L_mainlg,M,theta,c)
    S_noselg = Strouhal_number(f,L_noselg,M,theta,c)

    F_cleanwing = dimensionless_empirical_spectral_func(S_cleanwing,n_main,"Clean wing")
    F_slats = dimensionless_empirical_spectral_func(S_slats,n_main,"Slats")
    F_flaps = dimensionless_empirical_spectral_func(S_flaps,n_main,"Flaps")
    F_mainlg = dimensionless_empirical_spectral_func(S_mainlg,n_main,"Main landing gear")
    F_noselg = dimensionless_empirical_spectral_func(S_noselg,n_nose,"Main landing gear")

    D_cleanwing = directivity(theta,phi,df,"Clean wing")
    D_slats = directivity(theta,phi,df,"Slats")
    D_flaps = directivity(theta,phi,df,"Flaps")
    D_mainlg = directivity(theta,phi,df,"Main landing gear")
    D_noselg = directivity(theta,phi,df,"Main landing gear")

    pe_sq_cleanwing = pe_squared(rho_inf,c,P_cleanwing,D_cleanwing,F_cleanwing,r,M,theta)
    pe_sq_slats = pe_squared(rho_inf,c,P_slats,D_slats,F_slats,r,M,theta)
    pe_sq_flaps = pe_squared(rho_inf,c,P_flaps,D_flaps,F_flaps,r,M,theta)
    pe_sq_mainlg = pe_squared(rho_inf,c,P_mainlg,D_mainlg,F_mainlg,r,M,theta)
    pe_sq_noselg = pe_squared(rho_inf,c,P_noselg,D_noselg,F_noselg,r,M,theta)

    pe_sq[:,i] = (pe_sq_cleanwing + pe_sq_slats + pe_sq_flaps + 2*pe_sq_mainlg + 3*pe_sq_noselg)[0]

print(pe_sq)
print(pe0**2)
SPL = SPL(pe_sq,pe0)
#print(SPL)
dL_A = dL_A(f)

L_A = L_A(SPL,dL_A)

print(np.max(L_A))

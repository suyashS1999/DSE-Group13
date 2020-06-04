import numpy as np

def sweep_x(x,lmb_le,Cr,Ct,b_total):
	
	sweep_arb = np.arctan((np.tan(lmb_le)) - (Cr*x/(b_total*0.5)) + (Ct*x/(b_total*0.5)))
	
	return sweep_arb
	

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



def e(AR,sweep_le):
	
	e = (4.61*(1-(0.045*(AR**0.68)))*((np.cos(sweep_le))**0.15)) - 3.1
	
	return e


def CL_alpha_DATCOM(A,M,a0,sweep_half):
	
	k = a0/(2*np.pi)
	a = (2*np.pi*A)/(2+np.sqrt((((A**2)*(1-(M**2))/k**2)*(1+((np.tan(sweep_half)**2)/(1-(M**2)))))+4))
	
	return a

def CD0_wing(S_wet_ratio,t_c_avg,x_c_m,sweep_LE,C_r_m,C_t_m,C_r,b_total,rho_cruise,V_cruise,mu_cruise,k_wing,M_cruise,laminar_flow,turb_flow):
	

	sweep_max_thickness = sweep_x(x_c_m,sweep_LE,C_r_m, C_t_m, b_total)

	Re_cruise_max = Re_transonic(rho_cruise, V_cruise, C_r,mu_cruise, k_wing, M_cruise)


	Cf_total = laminar_flow*Cf_laminar(Re_cruise_max) + turb_flow*Cf_turb(M_cruise, Re_cruise_max)

	FF_wing = Form_factor_wing(t_c_avg, x_c_m, M_cruise, sweep_max_thickness)

	


	CD_0_wing = ((S_wet_ratio*Cf_total*FF_wing))
	
	return CD_0_wing
	
	

def delta_wave_drag(M_cruise,t_c_stream,sweep_LE,CL_cruise,M_crit_airfoil):
	
	M_dd_cruise = M_dd(t_c_stream, sweep_LE, CL_cruise)
	M_crit_wing = M_crit_airfoil/np.cos(sweep_LE)
	
	if M_crit_wing>M_cruise:
		delta_cd = 0
		print("Gooood")

	elif M_crit_wing<=M_cruise and M_cruise<=M_dd_cruise:
		delta_cd = 0.002*((1+(2.5*((M_dd_cruise-M_cruise)/0.05)))**-1)
		print("Okeii")
	elif M_dd_cruise<M_cruise:
		delta_cd = 0.002*((1+(((M_dd_cruise-M_cruise)/0.05)))**25)
		print("Baaaaaaddd")
	
	return delta_cd





# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 16:59:13 2020

@author: sanjay
"""


import numpy as np


# Cp_min = -0.43

gam = 1.4
M_crit = 0.72

cpm = -0.43


# def M_cr_calc(M_cr):
# 	
# 	gam = 1.4
# 	
#  	cpm = -0.31
# 	p1 = ((2+((gam-1)*(M_cr**2)))/(gam+1))**(gam/(gam-1))
# 	
# 	val = (cpm/np.sqrt(1-(M_cr**2))) -(2/(gam*(M_cr**2)))*(p1-1)
# 	
# 	return val



err = 1
count = 0


# while abs(err)>1e-3:
# 	print("Iteration ", count)
# 	val_1 = modif_cp(Cp_min, M_crit)
# 	val_2 = M_cr_calc(gam, M_crit)
# 	
# 	err = abs(val_1)-abs(val_2)
# 	if err>0 or err<0.1:
# 		M_crit = M_crit -0.01
# 		
#     elif err>0.1 or err<0.01:
# 		M_crit = M_crit - 0.001
# 		
# 	elif err<0 or err<0.1:
# 	
# 	else:
# 		M_crit = M_crit + 0.01
# 	
# 	count = count +1


def bisection(f,a,b,N):
    '''
    Parameters
    ----------
    f : function for which solution f(x)=0 is needed
    a,b : The interval in which to search for a solution. 
    N : number of iterations 

    '''
    if f(a)*f(b) >= 0:
        print("Bisection method fails.")
        return None
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = f(m_n)
        if f(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n
        else:
            print("Bisection method fails.")
            return None
    return (a_n + b_n)/2




"""Change cpm to particular airfoil value"""
def M_cr_calc(M_cr):
	
	gam = 1.4
	
	cpm = -0.43  # min cp for airfoil
	
	p1 = ((2+((gam-1)*(M_cr**2)))/(gam+1))**(gam/(gam-1))
	
	val = (cpm/np.sqrt(1-(M_cr**2))) -(2/(gam*(M_cr**2)))*(p1-1)
	
	return val



N = 100	
print("M_crit =", bisection(M_cr_calc,0.2,0.99,N))
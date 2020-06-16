# -*- coding: utf-8 -*-
"""

"""

import Aeroelas_analysis as aero
import sympy as sym
import numpy as np

x = sym.Symbol("x")
lmb = sym.Symbol("t")
V = 80
rho = 1.225
s= 7.5  #aero.half_span

mu = (aero.GJ* (s**2))/((s**2)*4*aero.EI)

a1 = aero.A_mat[0,0]*(lmb**2) + rho*aero.B_mat[0,0]*V*lmb + (rho*aero.C_mat[0,0] + x)*(V**2)

a2 = aero.A_mat[0,1]*(lmb**2) + rho*aero.B_mat[0,1]*V*lmb + (rho*aero.C_mat[0,1])*(V**2)

a3 = aero.A_mat[1,0]*(lmb**2) + rho*aero.B_mat[1,0]*V*lmb + (rho*aero.C_mat[1,0])*(V**2)

a4 = aero.A_mat[1,1]*(lmb**2) + rho*aero.B_mat[1,1]*V*lmb + (rho*aero.C_mat[1,1] + (mu*x))*(V**2)


B = np.array([[a1,a2],[a3,a4]])

B = sym.Matrix(B)
check = B.det()
check = check.expand()


b4 = check.coeff(lmb,n=4)
b3 = check.coeff(lmb,n=3)
b2 = check.coeff(lmb,n=2)#+ check.coeff(x).coeff(lmb**2)
b1 = check.coeff(lmb,n =1)
b0 = check.coeff(lmb,n=0)

result = (b4*(b1**2)) - (b1*b2*b3) + (b0*(b3**2))


sol = sym.solve(result)


# E = aero.E_mat

# tryy = E[0,0]/abs(sol[0])
# V_crit_1 = sym.sqrt(tryy)

# V_crit_2 = sym.sqrt(E[1,1]/(mu*np.absolute(sol[1])))

# V_crit_3 = sym.sqrt(E[0,0]/abs(sol[1]))

# V_crit_4 = sym.sqrt(E[1,1]/(mu*np.absolute(sol[0])))